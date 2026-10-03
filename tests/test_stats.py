# -*- coding: utf-8 -*-
"""src/stats.py 的单元测试,全部使用手工可验证的期望值。

运行方式(两种皆可):
    python tests/test_stats.py          # 独立运行,无需 pytest
    python -m pytest tests -q           # 若已安装 pytest

所有期望值均由公式手工推导,推导过程写在各测试的注释中。
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.stats import (  # noqa: E402
    ContingencyTable,
    bcpnn_ic,
    chi2_sf,
    contingency_table_from_counts,
    evaluate_pair,
    normal_cdf,
    normal_ppf,
    proportional_reporting_ratio,
    reporting_odds_ratio,
)

TOL = 1e-6


def approx(a: float, b: float, tol: float = TOL) -> bool:
    """相对容差比较;tol 作为绝对下限,避免 b 接近 0 时要求不可能的相对精度。"""
    return abs(a - b) <= max(tol, tol * abs(b))


# ---------------------------------------------------------------------------
# 基础分布函数
# ---------------------------------------------------------------------------
def test_normal_cdf_known_values():
    # Phi(0) = 0.5;Phi(1.96) ≈ 0.9750021;Phi(-1.96) ≈ 0.0249979
    assert approx(normal_cdf(0.0), 0.5, 1e-7)
    assert approx(normal_cdf(1.96), 0.9750021, 1e-5)
    assert approx(normal_cdf(-1.96), 0.0249979, 1e-5)
    assert approx(normal_cdf(1.0), 0.8413447, 1e-5)


def test_normal_ppf_known_values():
    # 分位数应为 normal_cdf 的数值逆
    assert approx(normal_ppf(0.975), 1.959964, 1e-5)
    # 中位数需要绝对精度:1e-9 已足够严格(近似式 + Halley 精修后偏差约 1e-9)
    assert abs(normal_ppf(0.5)) < 1e-8
    assert approx(normal_ppf(0.025), -1.959964, 1e-5)
    # 往返一致性(含尾部区域)
    for p in (1e-4, 0.01, 0.3, 0.9, 0.999):
        assert approx(normal_cdf(normal_ppf(p)), p, 1e-6)


def test_chi2_sf_known_values():
    # df=1 时,chi2=3.841459 对应上尾概率 0.05
    assert approx(chi2_sf(3.841459, 1), 0.05, 1e-5)
    assert approx(chi2_sf(6.634897, 1), 0.01, 1e-4)
    assert approx(chi2_sf(0.0, 1), 1.0, 1e-12)
    # df=2 时 chi2=5.991465 对应 0.05
    assert approx(chi2_sf(5.991465, 2), 0.05, 1e-4)
    # chi2=1, df=1:p = 2*(1-Phi(1)) = 0.3173105
    assert approx(chi2_sf(1.0, 1), 0.3173105, 1e-5)


# ---------------------------------------------------------------------------
# PRR
# ---------------------------------------------------------------------------
def test_prr_hand_computed():
    """手算: a=20, b=80, c=30, d=370, N=500。

    PRR = [20/100] / [30/400] = 0.2 / 0.075 = 8/3 = 2.6666667
    chi2 = N * (a*d - b*c)^2 / [ (a+b)(c+d)(a+c)(b+d) ]
         = 500 * (20*370 - 80*30)^2 / (100 * 400 * 50 * 450)
         = 500 * 5000^2 / 900000000 = 12500000000 / 900000000
         = 125/9 = 13.8888889
    (已用独立定义的 Pearson 卡方 sum (O-E)^2/E 交叉验证,两者一致。)
    SE(ln PRR) = sqrt(1/20 - 1/100 + 1/30 - 1/400)
               = sqrt(0.05 - 0.01 + 0.0333333 - 0.0025)
               = sqrt(0.0708333) = 0.2661458
    ln PRR = 0.9808293
    CI 低 = exp(0.9808293 - 0.5216458) = exp(0.4591835) = 1.5827826
    CI 高 = exp(0.9808293 + 0.5216458) = exp(1.5024751) = 4.4927909
    p 值 = P(chi2_1 > 13.8888889) = 0.00019394
    """
    t = ContingencyTable(20, 80, 30, 370)
    r = proportional_reporting_ratio(t)

    assert r["corrected"] is False
    assert approx(r["prr"], 8.0 / 3.0, 1e-9)
    assert approx(r["chi2"], 125.0 / 9.0, 1e-9)
    assert approx(r["ci_low"], 1.5827826, 1e-6)
    assert approx(r["ci_high"], 4.4927909, 1e-6)
    assert approx(r["p_value"], 0.00019394, 1e-7)


def test_prr_no_association():
    """无关联:比例相同,PRR 应恰为 1,chi2 为 0。"""
    t = ContingencyTable(10, 90, 10, 90)
    r = proportional_reporting_ratio(t)
    assert approx(r["prr"], 1.0, 1e-12)
    assert approx(r["chi2"], 0.0, 1e-12)
    assert approx(r["p_value"], 1.0, 1e-9)


def test_prr_zero_cell_uses_haldane():
    """a=0 时必须启用 +0.5 校正且不抛异常。

    校正后表为 (0.5, 80.5, 30.5, 370.5):
    PRR = [0.5/81] / [30.5/401] = 0.00617284 / 0.07605985 = 0.0811576
    """
    t = ContingencyTable(0, 80, 30, 370)
    r = proportional_reporting_ratio(t)
    assert r["corrected"] is True
    assert r["prr"] is not None and math.isfinite(r["prr"])
    assert approx(r["prr"], (0.5 / 81.0) / (30.5 / 401.0), 1e-9)

    # 校正不会被误用于非零表:表本身不变
    t2 = ContingencyTable(5, 5, 5, 5)
    assert t2.corrected() is t2


def test_prr_undefined_when_no_exposure():
    """药物完全无报告时统计量应为 None 而非崩溃。"""
    r = proportional_reporting_ratio(ContingencyTable(0, 0, 30, 370))
    assert r["prr"] is None, "分母为 0 时 PRR 必须为 None"


def test_prr_undefined_for_empty_table():
    """全零表不能返回 PRR=1 这一具有误导性的'无关联'结论。"""
    r = proportional_reporting_ratio(ContingencyTable(0, 0, 0, 0))
    assert r["prr"] is None
    assert r["chi2"] is None

    r2 = reporting_odds_ratio(ContingencyTable(0, 0, 0, 0))
    assert r2["ror"] is None


# ---------------------------------------------------------------------------
# ROR
# ---------------------------------------------------------------------------
def test_ror_hand_computed():
    """手算: a=20, b=80, c=30, d=370。

    ROR = (20*370) / (80*30) = 7400 / 2400 = 3.0833333
    SE(ln ROR) = sqrt(1/20 + 1/80 + 1/30 + 1/370)
               = sqrt(0.05 + 0.0125 + 0.0333333 + 0.0027027)
               = sqrt(0.0985360) = 0.3139051
    ln ROR = 1.1260113
    CI 低 = exp(1.1260113 - 0.6152540) = exp(0.5107573) = 1.6665420
    CI 高 = exp(1.1260113 + 0.6152540) = exp(1.7412653) = 5.7044000
    """
    t = ContingencyTable(20, 80, 30, 370)
    r = reporting_odds_ratio(t)

    assert r["corrected"] is False
    assert approx(r["ror"], 3.08333333, 1e-6)
    assert approx(r["se_ln_ror"], 0.3139051, 1e-5)
    assert approx(r["ci_low"], 1.6665420, 1e-4)
    assert approx(r["ci_high"], 5.7044000, 1e-4)


def test_ror_null_is_one():
    """交叉乘积相等时 ROR = 1。"""
    r = reporting_odds_ratio(ContingencyTable(10, 20, 20, 40))
    assert approx(r["ror"], 1.0, 1e-12)
    assert r["ci_low"] < 1.0 < r["ci_high"]


def test_ror_zero_cell_haldane():
    """b=0 时用 +0.5 校正:表变 (20.5, 0.5, 30.5, 370.5)。

    ROR = (20.5*370.5)/(0.5*30.5) = 7595.25 / 15.25 = 498.04918
    """
    r = reporting_odds_ratio(ContingencyTable(20, 0, 30, 370))
    assert r["corrected"] is True
    assert approx(r["ror"], (20.5 * 370.5) / (0.5 * 30.5), 1e-6)
    assert r["ci_low"] > 1.0


def test_ror_symmetric_inversion():
    """ROR 具有倒数对称性:交换行后 ROR 取倒数。"""
    a = reporting_odds_ratio(ContingencyTable(20, 80, 30, 370))["ror"]
    b = reporting_odds_ratio(ContingencyTable(30, 370, 20, 80))["ror"]
    assert approx(a * b, 1.0, 1e-9)


# ---------------------------------------------------------------------------
# BCPNN IC
# ---------------------------------------------------------------------------
def test_bcpnn_ic_hand_computed():
    """手算: a=20, b=80, c=30, d=370, N=500。

    E = (a+b)(a+c)/N = 100*50/500 = 10
    IC(近似观测量) = log2(20/10) = 1.0

    E(IC) = log2( (a+g11)(N+a0)(N+b0) / [ (N+g)(a+b+a1)(a+c+b1) ] )
          = log2( 20.5 * 501 * 501 / (501 * 100.5 * 50.5) )
          = log2( 5144200.5 / 2542700.25 ) = log2(2.02313...) = 1.0169556

    V(IC) = 1/(ln2)^2 * [ (N-a+g-g11)   / ((a+g11)(1+N+g))
                          + (N-(a+b)+a0-a1) / ((a+b+a1)(1+N+a0))
                          + (N-(a+c)+b0-b1) / ((a+c+b1)(1+N+b0)) ]
          = 2.0813681 * [ 480.5/(20.5*502) + 400.5/(100.5*502) + 450.5/(50.5*502) ]
          = 0.1506915
    sd = sqrt(0.1506915) = 0.3881900
    IC025 = 1.0169556 - 1.96*0.3881900 = 0.2561032
    IC975 = 1.0169556 + 1.96*0.3881900 = 1.7778081
    """
    t = ContingencyTable(20, 80, 30, 370)
    r = bcpnn_ic(t)

    assert approx(r["expected"], 10.0, 1e-9)
    assert approx(r["ic_raw"], 1.0, 1e-9)
    assert approx(r["ic"], 1.0169556, 1e-6)
    assert approx(r["ic_var"], 0.1506915, 1e-6)
    assert approx(r["ic025"], 0.2561032, 1e-6)
    assert approx(r["ic975"], 1.7778081, 1e-6)


def test_bcpnn_ic_zero_when_at_expectation():
    """观测等于期望时 IC 的近似值应为 0,先验收缩使其略偏正。

    取 a=10, b=90, c=90, d=810:N=1000,E = 100*100/1000 = 10 = a。
    ic_raw = log2(10/10) = 0
    先验 E(IC) = log2( 10.5*1001*1001 / (1001 * 100.5 * 90.5) ) = 0.0574403
    收缩量约 0.06,量级很小,这是弱信息先验在 a=10 时的正常表现。
    """
    t = ContingencyTable(10, 90, 90, 810)
    r = bcpnn_ic(t)
    assert approx(r["expected"], 10.0, 1e-9)
    assert approx(r["ic_raw"], 0.0, 1e-12)
    # 与近似值同号且量级很小(先验收缩不会改变"无关联"的结论)
    assert abs(r["ic"]) < 0.1


def test_bcpnn_ic025_ordering():
    """IC025 必须小于等于 IC,IC975 必须大于等于 IC。"""
    r = bcpnn_ic(ContingencyTable(20, 80, 30, 370))
    assert r["ic025"] < r["ic"] < r["ic975"]


def test_bcpnn_ic_undefined_without_cases():
    """a=0 时 IC 无定义(应为 None),不应抛异常。"""
    r = bcpnn_ic(ContingencyTable(0, 80, 30, 370))
    assert r["ic"] is None
    assert r["ic025"] is None


def test_bcpnn_ic_increases_with_signal():
    """不相称程度增强时 IC 与 IC025 应单调上升。"""
    weak = bcpnn_ic(ContingencyTable(12, 88, 30, 370))["ic025"]
    strong = bcpnn_ic(ContingencyTable(40, 60, 30, 370))["ic025"]
    assert strong > weak


# ---------------------------------------------------------------------------
# 列联表构造与信号判据
# ---------------------------------------------------------------------------
def test_contingency_table_from_counts():
    """由 a=20, n1=100, m1=50, N=500 反推 b=80, c=30, d=370。"""
    t = contingency_table_from_counts(20, 100, 50, 500)
    assert (t.a, t.b, t.c, t.d) == (20.0, 80.0, 30.0, 370.0)
    assert t.n == 500.0
    assert t.n_drug == 100.0
    assert t.n_event == 50.0


def test_contingency_table_negative_guard():
    """口径不一致导致的负格应被截断为 0,而不是静默产生错误值。"""
    t = contingency_table_from_counts(60, 100, 50, 500)  # m1=50 < a=60,不合理输入
    assert t.c == 0.0
    assert min(t.a, t.b, t.c, t.d) >= 0.0


def test_evaluate_pair_signal_flags():
    """a=20,b=80,c=30,d=370:PRR=2.667>=2、chi2=13.889>=4、a=20>=3,三条 EU 条件均成立;
    ROR 下限 1.667 > 1 成立;IC025 = 0.256 > 0 成立。故三判据同时命中。
    """
    t = ContingencyTable(20, 80, 30, 370)
    r = evaluate_pair("TESTDRUG", "TESTREACTION", t)

    assert r.is_signal_prr is True
    assert r.is_signal_ror is True
    assert r.is_signal_ic is True
    assert r.is_signal is True
    assert r.a == 20


def test_prr_criterion_requires_chi2():
    """PRR 判据必须同时满足 chi2 >= 4,仅有高 PRR 并不足够。

    取 a=3, b=23, c=1, d=59 (N=86):
      PRR = (3/26) / (1/60) = 0.1153846 / 0.0166667 = 6.9230769  (远大于 2)
      chi2 = 86*(3*59 - 23*1)^2 / (26 * 60 * 4 * 82)
           = 86 * 154^2 / 511680 = 86*23716/511680 = 2039576/511680
           = 3.9860381           (< 4,故 EU 判据不成立)
    样本量小时的这一情形正说明:单看 PRR 的大小会产生误判。
    """
    t = ContingencyTable(3, 23, 1, 59)
    r = evaluate_pair("D", "E", t)
    assert approx(r.prr, 6.9230769, 1e-6)
    assert r.chi2 < 4.0
    assert r.is_signal_prr is False, "chi2<4 时 PRR 判据不应成立"


def test_prr_criterion_requires_three_cases():
    """PRR 与 chi2 均达标时,a < 3 仍不构成 EU 判据。

    a=2, b=198, c=1, d=799:
      PRR = (2/200) / (1/800) = 0.01 / 0.00125 = 8.0  (>= 2)
      chi2 = 1000*(2*799 - 198*1)^2 / (200 * 800 * 3 * 997)
           = 1000*1400^2 / 478560000 = 1960000000/478560000 = 4.0956  (>= 4)
    两项数值条件都满足,唯一不满足的是最少 3 例的要求。
    """
    t = ContingencyTable(2, 198, 1, 799)
    r = evaluate_pair("D", "E", t)
    assert approx(r.prr, 8.0, 1e-9)
    assert r.chi2 >= 4.0
    assert r.is_signal_prr is False, "a=2 < 3,PRR 判据不应成立"


def test_evaluate_pair_strong_signal():
    """构造一个三判据同时成立的组合:a=50,b=50,c=50,d=850,N=1000。

    PRR = 0.5 / 0.0555556 = 9.0
    chi2 = 1000*(50*850-50*50)^2/(100*900*100*900) = 1000*1600000000/8100000000 = 197.53
    ROR = (50*850)/(50*50) = 17.0
    E = 100*100/1000 = 10,E(IC) ≈ log2((50.5*1001^2)/(1001*100.5*100.5)) ≈ 2.32
    """
    t = ContingencyTable(50, 50, 50, 850)
    r = evaluate_pair("STRONG", "RARE EVENT", t)

    assert approx(r.prr, 9.0, 1e-9)
    assert approx(r.ror, 17.0, 1e-9)
    assert r.chi2 > 4
    assert r.is_signal_prr and r.is_signal_ror and r.is_signal_ic
    assert r.is_signal is True


def test_evaluate_pair_n_below_three_blocks_prr():
    """见 test_prr_criterion_requires_three_cases,此处保留独立断言以防回归。"""
    t = ContingencyTable(2, 198, 1, 799)
    r = evaluate_pair("D", "E", t)
    assert approx(r.prr, 8.0, 1e-9)
    assert r.chi2 >= 4.0
    assert r.is_signal_prr is False, "a=2 < 3,PRR 判据不应成立"


def test_require_all_criteria():
    """require_all_criteria=True 时,只有三判据同时成立才标记。
    a=20,b=80,c=30,d=370 三判据均成立,故标记为真。"""
    t = ContingencyTable(20, 80, 30, 370)
    r = evaluate_pair("D", "E", t, require_all_criteria=True)
    assert r.is_signal is True

    # 仅部分判据成立的表:a=3,b=23,c=1,d=59 的 chi2 未达标
    t2 = ContingencyTable(3, 23, 1, 59)
    r2 = evaluate_pair("D", "E", t2, require_all_criteria=True)
    assert r2.is_signal is False


def test_degenerate_table_no_crash():
    """全零表不应引发异常,所有统计量应为 None 或有限值。"""
    t = ContingencyTable(0, 0, 0, 0)
    r = evaluate_pair("D", "E", t)
    assert r.is_signal is False
    assert r.prr is None and r.ror is None and r.ic is None


def test_all_results_finite():
    """对一组随机化的表格做健全性检查:不出现 nan / inf。"""
    import random

    rng = random.Random(20260101)
    for _ in range(300):
        t = ContingencyTable(
            rng.randint(0, 500), rng.randint(0, 500),
            rng.randint(0, 500), rng.randint(0, 5000),
        )
        r = evaluate_pair("D", "E", t)
        for name in ("prr", "chi2", "ror", "ic", "ic025"):
            v = getattr(r, name)
            if v is not None:
                assert math.isfinite(v), f"{name} 非有限值: {v} (表 {t})"


# ---------------------------------------------------------------------------
# 独立运行入口(无 pytest 时使用)
# ---------------------------------------------------------------------------
def _run_standalone() -> int:
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f)]
    passed, failed = 0, []
    for name, fn in tests:
        try:
            fn()
            passed += 1
            print(f"  PASS  {name}")
        except AssertionError as e:
            failed.append(name)
            print(f"  FAIL  {name}: {e}")
        except Exception as e:  # noqa: BLE001
            failed.append(name)
            print(f"  ERROR {name}: {type(e).__name__}: {e}")
    print()
    print(f"共 {len(tests)} 个测试,通过 {passed},失败 {len(failed)}")
    if failed:
        print("失败用例:" + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError, OSError):
        pass
    raise SystemExit(_run_standalone())
