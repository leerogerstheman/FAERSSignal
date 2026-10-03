"""Markdown 报告生成。

报告固定包含:方法学、2x2 列联表定义、信号判据、结果表(按信号强度排序)、
以及局限性说明。所有数字均直接来自本次运行的计算结果。
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Optional

from .stats import SignalResult, SIGNAL_CRITERIA_TEXT, sort_key_signal_strength

__all__ = ["render_markdown", "write_markdown"]


def _fmt(value: Optional[float], digits: int = 3) -> str:
    """数值格式化:None 统一显示为 "-",避免出现 nan/None 混排。"""
    if value is None:
        return "-"
    try:
        f = float(value)
    except (TypeError, ValueError):
        return "-"
    if f != f:  # NaN
        return "-"
    if f in (float("inf"), float("-inf")):
        return "inf" if f > 0 else "-inf"
    return f"{f:.{digits}f}"


def _fmt_ci(low: Optional[float], high: Optional[float], digits: int = 3) -> str:
    if low is None or high is None:
        return "-"
    return f"{_fmt(low, digits)}-{_fmt(high, digits)}"


def _result_row(r: SignalResult) -> str:
    flags = []
    if r.is_signal_prr:
        flags.append("PRR")
    if r.is_signal_ror:
        flags.append("ROR")
    if r.is_signal_ic:
        flags.append("IC")
    flag = "+".join(flags) if flags else "no"
    drug = r.drug.replace("|", "\\|")
    reaction = r.reaction.replace("|", "\\|")
    return (
        f"| {drug} | {reaction} | {int(r.a)} | {_fmt(r.prr)} | {_fmt_ci(r.prr_ci_low, r.prr_ci_high)} "
        f"| {_fmt(r.chi2, 2)} | {_fmt(r.ror)} | {_fmt_ci(r.ror_ci_low, r.ror_ci_high)} "
        f"| {_fmt(r.ic, 2)} | {_fmt(r.ic025, 2)} | {flag} |"
    )


def render_markdown(
    results: Iterable[SignalResult],
    drugs: Iterable[str],
    reactions: List[str],
    grand_total: int,
    top_n: int = 40,
    generated_at: Optional[str] = None,
    data_last_updated: Optional[str] = None,
    notes: Optional[List[str]] = None,
) -> str:
    """渲染完整 Markdown 报告。"""
    results = list(results)
    drugs = list(drugs)
    generated_at = generated_at or datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    signals = [r for r in results if r.is_signal]
    signals.sort(key=sort_key_signal_strength)
    strong = [r for r in signals if r.is_signal_prr and r.is_signal_ror and r.is_signal_ic]

    lines: List[str] = []
    A = lines.append

    A("# FAERS 药物警戒信号检测报告")
    A("")
    A(f"- 生成时间:{generated_at}")
    A(f"- 数据来源:openFDA `drug/event` API(FAERS 自发呈报数据)")
    if data_last_updated:
        A(f"- 数据源标注的最新更新时间:{data_last_updated}")
    A(f"- 分析药物({len(drugs)} 个):{', '.join(drugs)}")
    A(f"- 分析反应数:{len(reactions)} 个 MedDRA 首选术语(PT)")
    A(f"- 药物-反应组合总数:{len(results)};命中信号:{len(signals)};三条判据同时满足:{len(strong)}")
    A(f"- 全库药物-反应报告对总数 N = {grand_total:,}(见方法学中的口径说明)")
    A("")

    A("## 1. 方法学")
    A("")
    A("本研究采用**不相称性分析(disproportionality analysis)**,这是药物警戒中筛查")
    A("疑似安全信号的常规定量方法。其核心思想是:如果某药物与某不良事件之间不存在")
    A("特殊的关联,那么该事件在\"该药物报告\"中的占比应当与在其他药物报告中的占比接近;")
    A("若前者显著更高,则提示可能存在待评估的信号。")
    A("")
    A("对每一个\"药物-反应\"组合,按下列方式构造 2x2 列联表:")
    A("")
    A("```")
    A("                    目标反应 a        其他反应         合计")
    A("    目标药物           a               b            n1 = a+b")
    A("    其他药物           c               d            n0 = c+d")
    A("    合计            m1 = a+c        m0 = b+d         N = a+b+c+d")
    A("```")
    A("")
    A("其中 a/b/c/d 的含义为:")
    A("")
    A("- `a`:数据库中**该药物**同时报告了**该反应**的报告数")
    A("- `b`:该药物报告了其他反应的报告数")
    A("- `c`:其他药物报告了该反应的报告数")
    A("- `d`:其余报告数")
    A("- `N`:该数据库中\"药物-反应\"报告对的总数(即全库报告数的近似)")
    A("")
    A("计数字段来自 openFDA 的 `count=patient.reaction.reactionmeddrapt.exact` 聚合查询;")
    A("`a` 取自该药物反应分布中对应 PT 的计数,`n1 = a+b` 与 `m1 = a+c` 取自对应的")
    A("`meta.results.total`,因此不受单次聚合查询 `limit` 的截断影响。")
    A("")
    A("**连续性校正**:当 2x2 表中出现 0 单元格时,对数或比值会出现除零。")
    A("此时对四个单元格**同时加 0.5**(Haldane-Anscombe 校正),仅在这种情况下启用。")
    A("该做法使估计值向无效值(PRR=ROR=1)收缩,相关行在结果表中标注 `corrected`。")
    A("")

    A("## 2. 统计量定义")
    A("")
    A("### 2.1 PRR(比例报告比, Proportional Reporting Ratio)")
    A("")
    A("```")
    A("PRR = [ a / (a+b) ] / [ c / (c+d) ]")
    A("```")
    A("")
    A("标准误与 95% 置信区间(对数正态近似):")
    A("")
    A("```")
    A("SE(ln PRR) = sqrt( 1/a - 1/(a+b) + 1/c - 1/(c+d) )")
    A("95% CI     = exp( ln PRR +/- 1.96 * SE )")
    A("```")
    A("")
    A("卡方检验(2x2 表 Pearson 卡方,自由度 1):")
    A("")
    A("```")
    A("chi2 = N * (a*d - b*c)^2 / [ (a+b) * (c+d) * (a+c) * (b+d) ]")
    A("```")
    A("")
    A("### 2.2 ROR(报告比值比, Reporting Odds Ratio)")
    A("")
    A("```")
    A("ROR = (a / b) / (c / d) = a*d / (b*c)")
    A("SE(ln ROR) = sqrt( 1/a + 1/b + 1/c + 1/d )")
    A("95% CI     = exp( ln ROR +/- 1.96 * SE )")
    A("```")
    A("")
    A("### 2.3 BCPNN 的 IC(信息成分, Information Component)")
    A("")
    A("BCPNN(Bayesian Confidence Propagation Neural Network)由 WHO 乌普萨拉监测中心")
    A("(UMC)用于其信号检测。IC 定义为观测值与期望值之比的以 2 为底的对数:")
    A("")
    A("```")
    A("IC      = log2( a / E ),  E = (a+b) * (a+c) / N")
    A("E(IC)   = log2( (a + g11)(N + a0)(N + b0) / [ (N + g)(a + b + a1)(a + c + b1) ] )")
    A("V(IC)   = 1/(ln2)^2 * [ (N - a + g - g11) / ((a + g11)(1 + N + g))")
    A("                       + (N - (a+b) + a0 - a1) / ((a+b+a1)(1 + N + a0))")
    A("                       + (N - (a+c) + b0 - b1) / ((a+c+b1)(1 + N + b0)) ]")
    A("IC025   = E(IC) - 1.96 * sqrt(V(IC))")
    A("```")
    A("")
    A("先验超参数取 UMC 惯例的一组弱信息先验:g11=0.5, g=1.0, a0=1.0, b0=1.0, a1=0.5, b1=0.5。")
    A("在 N 较大时该式退化为 IC = log2(a/E)。**IC025 > 0** 表示在贝叶斯框架下")
    A("有 95% 的可信度认为观测数高于期望数,是常用的信号判据。")
    A("")

    A("## 3. 信号判据")
    A("")
    A(f"本报告采用的判据:{SIGNAL_CRITERIA_TEXT}。")
    A("")
    A("| 判据 | 阈值 | 出处 |")
    A("|---|---|---|")
    A("| PRR | PRR >= 2 且 chi2 >= 4 且 a >= 3 | 欧盟(EMA)对 PRR 的常规信号阈值 |")
    A("| ROR | ROR 的 95% CI 下限 > 1 | 荷兰 Lareb / 常规 ROR 判据 |")
    A("| IC | IC025 > 0 | WHO UMC BCPNN 判据 |")
    A("")
    A("结果表中 `判据` 列以 `PRR`/`ROR`/`IC` 标识各自独立的命中情况,")
    A("`PRR+ROR+IC` 表示三条判据同时满足(证据相对更强的组合),`no` 表示未命中任何判据。")
    A("")

    A("## 4. 结果")
    A("")
    if signals:
        A(f"共检出 **{len(signals)}** 个信号(任一判据命中),下表按信号强度(P RR 降序,"
          f"其次 IC025)列出前 {min(top_n, len(signals))} 个。")
        A("")
        A("| 药物 | 反应 (MedDRA PT) | a | PRR | PRR 95% CI | chi2 | ROR | ROR 95% CI | IC | IC025 | 判据 |")
        A("|---|---|---:|---:|---|---:|---:|---|---:|---:|---|")
        for r in signals[:top_n]:
            A(_result_row(r))
    else:
        A("本次分析未检出任何满足判据的信号。")
    A("")

    if strong:
        best = strong[0]
        A("### 4.1 三条判据同时命中的组合")
        A("")
        A("| 药物 | 反应 | a | PRR | ROR (95% CI) | IC025 |")
        A("|---|---|---:|---:|---|---:|")
        for r in strong[:15]:
            A(f"| {r.drug} | {r.reaction} | {int(r.a)} | {_fmt(r.prr)} "
              f"| {_fmt(r.ror)} ({_fmt_ci(r.ror_ci_low, r.ror_ci_high)}) | {_fmt(r.ic025, 2)} |")
        A("")
        A(f"其中不相称程度最高的是 **{best.drug} — {best.reaction}**:")
        A(f"PRR = {_fmt(best.prr)}(95% CI {_fmt_ci(best.prr_ci_low, best.prr_ci_high)})、"
          f"ROR = {_fmt(best.ror)}(95% CI {_fmt_ci(best.ror_ci_low, best.ror_ci_high)})、"
          f"IC025 = {_fmt(best.ic025, 2)},对应报告数 a = {int(best.a)}。")
        A("")

    A("### 4.2 全部组合明细")
    A("")
    all_sorted = sorted(results, key=sort_key_signal_strength)
    A("| 药物 | 反应 | a | b | c | d | PRR | ROR | IC | IC025 | 信号 |")
    A("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in all_sorted:
        A(f"| {r.drug} | {r.reaction} | {int(r.a)} | {int(r.b)} | {int(r.c)} | {int(r.d)} "
          f"| {_fmt(r.prr)} | {_fmt(r.ror)} | {_fmt(r.ic, 2)} | {_fmt(r.ic025, 2)} "
          f"| {'是' if r.is_signal else '否'} |")
    A("")

    if notes:
        A("### 4.3 运行日志与数据口径提示")
        A("")
        for n in notes:
            A(f"- {n}")
        A("")

    A("## 5. 局限性(Limitations)")
    A("")
    A("本节是解读上述结果的必要前提,任何结论都不得脱离这些限制:")
    A("")
    A("1. **自发呈报数据没有分母**。FAERS 收集的是自愿提交的报告,不存在暴露人群数")
    A("   (接受该药物治疗的患者总数),因此 PRR / ROR / IC **不是发病率、不是风险比、")
    A("   不是绝对风险**,也不能用来比较不同药物的安全性高低。")
    A("2. **不相称不等于因果**。本研究只能说明\"报告中出现了超出背景水平的聚集\",")
    A("   不能证明药物导致了该不良事件。信号需要经流行病学研究和临床评估后")
    A("   才能转化为安全性结论。本报告不构成任何医疗建议。")
    A("3. **报告偏倚(reporting bias)**。报告行为受事件严重程度、上市时间、媒体关注、")
    A("   监管行动、诉讼与文献报道等因素影响,新药与备受关注的药物更容易被报告。")
    A("4. **重复报告**。同一病例可能被患者、医生、企业分别提交,FAERS 中未完全去重,")
    A("   会同时抬高 a 与 b。")
    A("5. **适应症混杂(confounding by indication)**。某药物被用于的疾病本身可能")
    A("   导致该不良事件,从而产生虚假的不相称。")
    A("6. **事件与药物的分组口径**。反应使用 MedDRA 首选术语(PT)层级的精确匹配,")
    A("   同一医学概念的不同 PT(SMQs 未合并)会被拆散;药物名按 openFDA 规范化的")
    A("   `medicinalproduct` 精确匹配,复方制剂、商品名与拼写差异可能造成漏计或串计。")
    A("7. **计数口径为下限估计**。全库报告对总数 N 由按药物名聚合的前若干项求和得到")
    A("   (聚合接口 `limit` 上限为 1000),未覆盖报告数极少的尾部药物名;")
    A("   反应总报告数 m1 取自 `meta.results.total`,不受该限制。因此 N 略偏低,")
    A("   会使 IC 与卡方被轻微高估,PRR/ROR 的方向不受影响。")
    A("8. **多次比较**。对成百上千个组合同时检验会产生大量假阳性,")
    A("   本报告未做多重比较校正(药物警戒筛查阶段的惯例)。")
    A("9. **统计方法本身是筛查工具**。PRR/ROR/IC 的设计目的是\"减少需要人工审阅的报告量\",")
    A("   而不是\"判定因果\"。阈值是人为约定的经验值,不同机构使用的阈值并不统一。")
    A("10. **数据时点**。openFDA 数据为滚动更新,复现时结果会随后续报告累积而略有变化;")
    A("    本报告的生成时间与数据源标注时间见文首。")
    A("")

    A("## 6. 参考文献")
    A("")
    A("- Evans SJW, Waller PC, Davis S. Use of proportional reporting ratios (PRRs) for")
    A("  signal generation from spontaneous adverse drug reaction reports.")
    A("  *Pharmacoepidemiol Drug Saf*. 2001;10(6):483-486.")
    A("- Rothman KJ, Lanes S, Sacks ST. The reporting odds ratio and its advantages over")
    A("  the proportional reporting ratio. *Pharmacoepidemiol Drug Saf*. 2004;13(8):519-523.")
    A("- Bate A, Lindquist M, Edwards IR, et al. A Bayesian neural network method for")
    A("  adverse drug reaction signal generation. *Eur J Clin Pharmacol*. 1998;54(4):315-321.")
    A("- Norén GN, Hopstadius J, Bate A. Shrinkage observed-to-expected ratios for robust")
    A("  and transparent large-scale pattern discovery. *Stat Methods Med Res*. 2013;22(1):57-69.")
    A("- 国家药品监督管理局药品评价中心.《药物警戒质量管理规范》相关技术要求。")
    A("- openFDA drug/event API 文档:https://open.fda.gov/apis/drug/event/")
    A("")
    A("---")
    A("")
    A("本报告由 FAERSSignal 自动生成,数据来源于美国 FDA openFDA 公开接口。")
    A("openFDA 免责声明:不得依据 openFDA 数据做出医疗决策;数据未经完整验证。")
    A("")

    return "\n".join(lines)


def write_markdown(content: str, out_path: str | Path) -> Path:
    """写出 UTF-8 Markdown 文件(不带 BOM)。"""
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(content, encoding="utf-8", newline="\n")
    return p
