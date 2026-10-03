"""命令行入口:对指定药物集合执行不相称性分析,输出 SQLite + CSV + Markdown。

典型用法::

    python faers_signal.py --drugs ASPIRIN,IBUPROFEN --top-reactions 30 \\
        --db data/faers.sqlite --csv output/signals.csv --md output/report.md

请求预算:每分析一个药物需要 2 次请求(反应分布 + 药物报告总数),再外加
每个反应 1 次请求(反应总报告数,结果在药物间复用)以及 1 次全库总数请求。
按默认 1 秒/请求的礼貌限流,3 个药物 x 30 个反应约需 95 秒。
"""

from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

from .openfda import OpenFDAClient, OpenFDAError
from .report import render_markdown, write_markdown
from .stats import (
    ContingencyTable,
    SignalResult,
    contingency_table_from_counts,
    evaluate_pair,
    sort_key_signal_strength,
)
from .store import SignalStore

__all__ = ["build_parser", "run_analysis", "main"]

# 默认药物集:覆盖不同药理类别、不同上市年代,便于展示信号的差异性
DEFAULT_DRUGS = ["ASPIRIN", "IBUPROFEN", "METFORMIN", "ATORVASTATIN", "SERTRALINE"]

# 需要过滤掉的非特异性/程序性反应术语:
# 这类术语反映的是报告流程或给药问题,而非医学意义上的不良事件,
# 若纳入会产生大量无临床意义的"强信号"。
EXCLUDED_REACTIONS = {
    "OFF LABEL USE", "PRODUCT USE IN UNAPPROVED INDICATION", "DRUG INEFFECTIVE",
    "PRODUCT QUALITY ISSUE", "PRODUCT SUBSTITUTION ISSUE", "INCORRECT DOSE ADMINISTERED",
    "DRUG ADMINISTRATION ERROR", "PRODUCT ADMINISTRATION ERROR", "WRONG TECHNIQUE IN DRUG DELIVERY PROCESS",
    "PRODUCT LABEL CONFUSION", "PRODUCT LABEL ISSUE", "EXPIRED PRODUCT ADMINISTERED",
    "PRODUCT MEASURED POTENCY ISSUE", "THERAPEUTIC RESPONSE UNEXPECTED",
    "INTENTIONAL PRODUCT USE ISSUE", "MEDICATION ERROR", "POOR QUALITY PRODUCT",
    "ADVERSE EVENT", "NO ADVERSE EVENT", "PRODUCT PACKAGING ISSUE",
}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="faers_signal",
        description="基于 openFDA 数据的药物警戒不相称性分析(PRR / ROR / BCPNN-IC)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--drugs", default=",".join(DEFAULT_DRUGS),
                   help="逗号分隔的药物名(openFDA 规范名),默认:" + ",".join(DEFAULT_DRUGS))
    p.add_argument("--top-reactions", type=int, default=30,
                   help="每个药物取报告数最高的前 N 个反应作为候选,默认 30")
    p.add_argument("--db", default="data/faers.sqlite", help="SQLite 输出路径")
    p.add_argument("--csv", default="output/signals.csv", help="CSV 输出路径")
    p.add_argument("--md", default="output/report.md", help="Markdown 报告输出路径")
    p.add_argument("--api-key", default=None, help="可选的 openFDA API key(提高限流上限)")
    p.add_argument("--sleep", type=float, default=1.0, help="请求间隔秒数,默认 1.0")
    p.add_argument("--grand-total", type=int, default=0,
                   help="直接指定全库报告对总数 N,>0 时跳过该次查询")
    p.add_argument("--timeout", type=float, default=45.0, help="单次请求超时秒数")
    p.add_argument("--top-n", type=int, default=40, help="报告中列出的信号条数上限")
    p.add_argument("--quiet", action="store_true", help="不打印进度")
    return p


def _fmt(v: Optional[float], d: int = 3) -> str:
    if v is None:
        return ""
    try:
        return f"{float(v):.{d}f}"
    except (TypeError, ValueError):
        return ""


def run_analysis(args: argparse.Namespace) -> int:
    """执行一次完整分析。返回 0 表示成功。"""
    drugs = [d.strip().upper() for d in args.drugs.split(",") if d.strip()]
    if not drugs:
        print("错误:未指定任何药物", file=sys.stderr)
        return 2

    verbose = not args.quiet
    client = OpenFDAClient(api_key=args.api_key, min_interval=args.sleep,
                           timeout=args.timeout, verbose=verbose)
    notes: List[str] = []
    started = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    print(f"=== FAERSSignal 分析开始 ({started}) ===")
    print(f"药物: {', '.join(drugs)}")
    print(f"每药取前 {args.top_reactions} 个反应作为候选")

    # ---- 全库报告对总数 N ----
    if args.grand_total > 0:
        grand_total = int(args.grand_total)
        notes.append(f"全库报告对总数 N 由参数指定:{grand_total:,}")
    else:
        print("[1/3] 获取全库药物-反应报告对总数 N ...")
        grand_total = client.grand_total()
        used = client.grand_total_limit
        notes.append(
            f"全库报告对总数 N = {grand_total:,},由按药物名聚合的前 {used} 项求和得到;"
            "openFDA 对该聚合端点拒绝 limit=1000(HTTP 403),故取可用的最大 limit。"
            "未覆盖报告数极少的尾部药物名,因此 N 为下限估计,会使 IC 与卡方被轻微高估。"
        )
    if grand_total <= 0:
        print("错误:无法获取全库总数,分析中止", file=sys.stderr)
        return 1
    print(f"      N = {grand_total:,}")

    # ---- 逐药物拉取反应分布 ----
    print(f"[2/3] 拉取各药物的反应分布(每药 1 次聚合查询)...")
    per_drug: Dict[str, List] = {}
    for drug in drugs:
        try:
            rows = client.count_reactions_for_drug(drug, limit=200)
        except OpenFDAError as e:
            print(f"      [跳过] {drug}: {e}")
            notes.append(f"药物 {drug} 查询失败,已跳过:{e}")
            continue
        if not rows:
            print(f"      [跳过] {drug}: openFDA 未返回任何报告")
            notes.append(f"药物 {drug} 无匹配报告,已跳过。")
            continue
        per_drug[drug] = rows
        print(f"      {drug}: {len(rows)} 个反应术语,最高 {rows[0].term} ({rows[0].count:,})")

    if not per_drug:
        print("错误:没有任何药物取得数据", file=sys.stderr)
        return 1

    # 候选反应:各药物报告数最高的前 N 个,按药物内排名交错合并,再去重
    candidates: List[str] = []
    seen = set()
    for rank in range(args.top_reactions):
        for drug in per_drug:
            rows = per_drug[drug]
            if rank < len(rows):
                term = rows[rank].term.upper()
                if term in EXCLUDED_REACTIONS or term in seen:
                    continue
                seen.add(term)
                candidates.append(term)
    if not candidates:
        print("错误:候选反应集合为空", file=sys.stderr)
        return 1

    # ---- 反应总报告数 m1(跨药物复用) ----
    print(f"[3/3] 查询 {len(candidates)} 个反应的全库报告数 m1 ...")
    event_totals: Dict[str, int] = {}
    for term in candidates:
        try:
            event_totals[term] = client.reaction_report_total(term)
        except OpenFDAError as e:
            print(f"      [跳过] {term}: {e}")
            event_totals[term] = 0

    # ---- 计算 ----
    results: List[SignalResult] = []
    for drug, rows in per_drug.items():
        try:
            drug_total = client.drug_report_total(drug)
        except OpenFDAError:
            drug_total = sum(r.count for r in rows)
        if drug_total <= 0:
            continue
        counts = {r.term.upper(): r.count for r in rows}
        for term in candidates:
            m1 = event_totals.get(term, 0)
            if m1 <= 0:
                continue
            a = counts.get(term, 0)
            table = contingency_table_from_counts(a, drug_total, m1, grand_total)
            results.append(evaluate_pair(drug, term, table))

    if not results:
        print("错误:未产生任何可计算的结果", file=sys.stderr)
        return 1

    results.sort(key=sort_key_signal_strength)
    n_signal = sum(1 for r in results if r.is_signal)
    n_strong = sum(1 for r in results
                   if r.is_signal_prr and r.is_signal_ror and r.is_signal_ic)
    notes.append(
        f"共计算 {len(results)} 个药物-反应组合,检出 {n_signal} 个信号"
        f"(其中 {n_strong} 个同时满足 PRR/ROR/IC 三条判据)。"
    )
    notes.append(
        f"已过滤 {len(EXCLUDED_REACTIONS)} 个非特异性/程序性术语"
        "(如 DRUG INEFFECTIVE、OFF LABEL USE、PRODUCT QUALITY ISSUE),"
        "这类术语反映报告流程或给药问题而非医学不良事件。"
    )
    notes.append(f"本次运行共发起 {client.request_count} 次 openFDA 请求,"
                 f"请求间隔 {args.sleep:.1f}s。")

    # ---- 输出 ----
    db_path = Path(args.db)
    with SignalStore(db_path) as store:
        run_id = store.start_run(drugs=list(per_drug), n_reactions=len(candidates),
                                 grand_total=grand_total, notes=" | ".join(notes))
        n_saved = store.save_results(run_id, results)
        store.finish_run(run_id)
    print(f"[输出] SQLite: {db_path} (run_id={run_id}, {n_saved} 行)")

    csv_path = Path(args.csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "drug", "reaction", "a", "b", "c", "d", "n_total", "expected",
        "prr", "prr_ci_low", "prr_ci_high", "chi2", "chi2_p",
        "ror", "ror_ci_low", "ror_ci_high", "ic", "ic025", "ic975", "ic_var",
        "corrected", "is_signal_prr", "is_signal_ror", "is_signal_ic", "is_signal",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in results:
            row = r.as_dict()
            for k in ("prr", "prr_ci_low", "prr_ci_high", "chi2", "chi2_p", "ror",
                      "ror_ci_low", "ror_ci_high", "ic", "ic025", "ic975", "ic_var",
                      "expected"):
                row[k] = "" if row.get(k) is None else round(float(row[k]), 6)
            w.writerow(row)
    print(f"[输出] CSV:    {csv_path}")

    md_content = render_markdown(
        results=results, drugs=list(per_drug), reactions=candidates,
        grand_total=grand_total, top_n=args.top_n,
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        notes=notes,
    )
    md_path = write_markdown(md_content, args.md)
    print(f"[输出] Markdown: {md_path}")

    print()
    print(f"=== 完成:{len(results)} 个组合,{n_signal} 个信号(强信号 {n_strong} 个)===")
    top = [r for r in results if r.is_signal][:5]
    for r in top:
        print(f"  {r.drug:14s} {r.reaction[:34]:34s} a={int(r.a):6d} "
              f"PRR={_fmt(r.prr)} ROR={_fmt(r.ror)} IC025={_fmt(r.ic025, 2)}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    # Windows 控制台默认 GBK,中文输出会抛 UnicodeEncodeError。
    # 尽量切到 UTF-8;失败时退化为 replace,绝不让编码问题中断分析。
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError, OSError):
        pass

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run_analysis(args)
    except KeyboardInterrupt:
        print("\n已中断", file=sys.stderr)
        return 130
    except OpenFDAError as e:
        print(f"openFDA 访问失败:{e}", file=sys.stderr)
        return 1
