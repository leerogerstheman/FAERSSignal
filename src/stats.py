"""不相称性分析(disproportionality analysis)统计量实现。

本模块不依赖 scipy / statsmodels,所有分布函数与区间估计均自实现。

术语与记号
----------
对某个"药物-反应"组合,基于自发呈报数据库的 2x2 列联表:

                目标反应 a          其他反应         合计
    目标药物      a                   b               n1 = a+b
    其他药物      c                   d               n0 = c+d
    合计          m1 = a+c            m0 = b+d        N  = a+b+c+d

其中:
    a = 目标药物的目标反应报告数
    b = 目标药物的其他反应报告数
    c = 其他药物的目标反应报告数
    d = 其他药物的其他反应报告数
    N = 数据库中全部"药物-反应"报告对的总数

注意事项:自发呈报数据无分母(暴露人群),因此以上统计量衡量的是
*报告不相称性(reporting disproportionality)*,不等同于风险比、发病率或因果关系。

连续性校正
----------
当 2x2 表中出现 0 单元格时,PRR 的比值、ROR 的对数以及 IC 的对数都会出现除零或
log(0)。本模块采用 Haldane-Anscombe 校正:对**全部四个单元格同时加 0.5**
(a,b,c,d → a+0.5, b+0.5, c+0.5, d+0.5),仅在原始表中存在 0 时启用。
该做法保持表格边缘合计的对称性,是流行病学中处理零单元格的常规手段,
但会使估计值偏向无效值(PRR/ROR 趋向 1),解读时需注意。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Iterable, Optional

__all__ = [
    "ContingencyTable",
    "SignalResult",
    "normal_cdf",
    "normal_ppf",
    "chi2_sf",
    "proportional_reporting_ratio",
    "reporting_odds_ratio",
    "bcpnn_ic",
    "contingency_table_from_counts",
    "evaluate_pair",
    "EU_SIGNAL_CRITERIA",
    "SIGNAL_CRITERIA_TEXT",
]

# 标准正态分布相关常量,避免依赖 scipy
_SQRT_2 = math.sqrt(2.0)
_SQRT_2PI = math.sqrt(2.0 * math.pi)

# EU(欧洲药品管理局)对 PRR 的信号判定阈值
EU_PRR_MIN = 2.0
EU_CHI2_MIN = 4.0
EU_N_MIN = 3

# 结构化描述,供报告与调用方复用,避免阈值散落在多处
EU_SIGNAL_CRITERIA = {
    "prr_min": EU_PRR_MIN,
    "chi2_min": EU_CHI2_MIN,
    "min_cases": EU_N_MIN,
    "ror_ci_low_min": 1.0,
    "ic025_min": 0.0,
}

SIGNAL_CRITERIA_TEXT = (
    "PRR >= 2.0 且 chi2 >= 4.0 且 a >= 3 (EU 标准);"
    "ROR 95% CI 下限 > 1.0;IC025 > 0"
)


# --------------------------------------------------------------------------
# 基础分布函数(自实现,无 scipy 依赖)
# --------------------------------------------------------------------------
def _erf(x: float) -> float:
    """误差函数,Abramowitz & Stegun 7.1.26 近似,绝对误差 < 1.5e-7。"""
    sign = 1.0 if x >= 0 else -1.0
    x = abs(x)
    t = 1.0 / (1.0 + 0.3275911 * x)
    y = 1.0 - (
        ((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
        + 0.254829592
    ) * t * math.exp(-x * x)
    return sign * y


def normal_cdf(x: float) -> float:
    """标准正态累积分布函数 Phi(x)。"""
    return 0.5 * (1.0 + _erf(x / _SQRT_2))


def normal_ppf(p: float) -> float:
    """标准正态分位数函数(逆 CDF),Acklam 有理逼近,相对误差 < 1.15e-9。"""
    if not (0.0 < p < 1.0):
        raise ValueError("normal_ppf 的 p 必须位于开区间 (0, 1)")

    a = (
        -3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02,
        1.383577518672690e02, -3.066479806614716e01, 2.506628277459239e00,
    )
    b = (
        -5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
        6.680131188771972e01, -1.328068155288572e01,
    )
    c = (
        -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
        -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00,
    )
    d = (
        7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
        3.754408661907416e00,
    )
    p_low, p_high = 0.02425, 1.0 - 0.02425

    if p < p_low:
        q = math.sqrt(-2.0 * math.log(p))
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )
    elif p <= p_high:
        q = p - 0.5
        r = q * q
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (
            ((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0
        )
    else:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / (
            (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0
        )

    # 一步 Halley 精修,把近似误差压到机器精度量级
    e = normal_cdf(x) - p
    u = e * _SQRT_2PI * math.exp(x * x / 2.0)
    x = x - u / (1.0 + x * u / 2.0)
    return x


def _lower_incomplete_gamma(s: float, x: float) -> float:
    """正则化下不完全伽马函数 P(s, x) = gamma(s, x) / Gamma(s)。

    级数展开用于 x < s+1,连分式展开用于其余情况(Numerical Recipes 6.2)。
    """
    if x <= 0.0:
        return 0.0
    gln = math.lgamma(s)
    eps = 3.0e-12
    itmax = 500

    if x < s + 1.0:
        ap = s
        total = 1.0 / s
        delta = total
        for _ in range(itmax):
            ap += 1.0
            delta *= x / ap
            total += delta
            if abs(delta) < abs(total) * eps:
                break
        return total * math.exp(-x + s * math.log(x) - gln)

    # 连分式:计算 Q(s, x) 再取补
    tiny = 1.0e-300
    b = x + 1.0 - s
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, itmax + 1):
        an = -i * (i - s)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    q = math.exp(-x + s * math.log(x) - gln) * h
    return 1.0 - q


def chi2_sf(chi2: float, df: int = 1) -> float:
    """卡方分布上尾概率 P(X > chi2),用于计算 p 值。

    df=1 时退化为 2*(1 - Phi(sqrt(chi2))),但统一走不完全伽马以支持任意自由度。
    """
    if chi2 <= 0.0:
        return 1.0
    if df <= 0:
        raise ValueError("自由度 df 必须为正整数")
    return 1.0 - _lower_incomplete_gamma(df / 2.0, chi2 / 2.0)


# --------------------------------------------------------------------------
# 2x2 列联表
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class ContingencyTable:
    """药物-反应 2x2 列联表的原始计数。"""

    a: float  # 目标药物 + 目标反应
    b: float  # 目标药物 + 其他反应
    c: float  # 其他药物 + 目标反应
    d: float  # 其他药物 + 其他反应

    def __post_init__(self) -> None:
        for name in ("a", "b", "c", "d"):
            v = getattr(self, name)
            if v < 0:
                raise ValueError(f"列联表单元格 {name} 不能为负: {v}")

    @property
    def n(self) -> float:
        """报告对总数 N。"""
        return self.a + self.b + self.c + self.d

    @property
    def n_drug(self) -> float:
        """目标药物的报告总数。"""
        return self.a + self.b

    @property
    def n_event(self) -> float:
        """目标反应的全部报告数。"""
        return self.a + self.c

    @property
    def has_zero_cell(self) -> bool:
        """是否存在零单元格(决定是否需要连续性校正)。"""
        return min(self.a, self.b, self.c, self.d) == 0

    def corrected(self) -> "ContingencyTable":
        """Haldane-Anscombe 校正:存在零单元格时四个格同时加 0.5。"""
        if not self.has_zero_cell:
            return self
        return ContingencyTable(self.a + 0.5, self.b + 0.5, self.c + 0.5, self.d + 0.5)

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class SignalResult:
    """单个药物-反应组合的全部分析结果。"""

    drug: str
    reaction: str
    a: float
    b: float
    c: float
    d: float
    n_total: float
    prr: Optional[float] = None
    prr_ci_low: Optional[float] = None
    prr_ci_high: Optional[float] = None
    chi2: Optional[float] = None
    chi2_p: Optional[float] = None
    ror: Optional[float] = None
    ror_ci_low: Optional[float] = None
    ror_ci_high: Optional[float] = None
    ic: Optional[float] = None
    ic025: Optional[float] = None
    ic975: Optional[float] = None
    ic_var: Optional[float] = None
    expected: Optional[float] = None
    corrected: bool = False
    is_signal_prr: bool = False
    is_signal_ror: bool = False
    is_signal_ic: bool = False

    @property
    def is_signal(self) -> bool:
        """三种判据任一成立即视为信号(EU PRR 判据为主判据)。"""
        return self.is_signal_prr or self.is_signal_ror or self.is_signal_ic

    def as_dict(self) -> dict:
        d = asdict(self)
        d["is_signal"] = self.is_signal
        return d


# --------------------------------------------------------------------------
# 统计量
# --------------------------------------------------------------------------
def proportional_reporting_ratio(table: ContingencyTable) -> dict:
    """PRR(比例报告比)及 chi2、95% 置信区间。

    公式:
        PRR = [a / (a+b)] / [c / (c+d)]

        chi2 = N * (a*d - b*c)^2 / [(a+b)(c+d)(a+c)(b+d)]   (无校正 Pearson,df=1)

        SE(ln PRR) = sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))
        95% CI     = exp(ln PRR +/- 1.96 * SE)

    返回 dict,含 prr / chi2 / p_value / ci_low / ci_high / corrected / n。
    出现零单元格时使用 Haldane-Anscombe(+0.5)校正后的表计算。
    """
    corrected = table.has_zero_cell
    t = table.corrected()

    out = {
        "prr": None,
        "chi2": None,
        "p_value": None,
        "ci_low": None,
        "ci_high": None,
        "corrected": corrected,
        "n": table.n,
    }

    # 分母为零(该药物或该反应完全没有报告)时统计量无定义。
    # 空表(四个格全为 0)也必须短路:否则校正后变成 0.5/0.5 的平凡比,
    # 会得到 PRR=1 这种"看似无关联"的伪结果。
    if table.n <= 0 or table.a + table.b <= 0 or (t.c + t.d) <= 0 or t.n_event <= 0:
        return out

    exp_drug = t.a / t.n_drug          # 目标药物中该反应的报告比例
    exp_other = t.c / (t.c + t.d)      # 其他药物中该反应的报告比例
    if exp_other <= 0:
        return out

    prr = exp_drug / exp_other
    out["prr"] = prr

    # Pearson 卡方(用未校正原表;若任一边际为 0 则无法计算)
    denom = table.n_drug * (table.c + table.d) * table.n_event * (table.b + table.d)
    if denom > 0 and table.n > 0:
        chi2 = table.n * (table.a * table.d - table.b * table.c) ** 2 / denom
        out["chi2"] = chi2
        out["p_value"] = chi2_sf(chi2, 1)

    # 对数正态近似置信区间
    try:
        se = math.sqrt(
            1.0 / t.a - 1.0 / t.n_drug + 1.0 / t.c - 1.0 / (t.c + t.d)
        )
        if prr > 0 and se > 0 and math.isfinite(se):
            ln = math.log(prr)
            out["ci_low"] = math.exp(ln - 1.96 * se)
            out["ci_high"] = math.exp(ln + 1.96 * se)
    except (ValueError, ZeroDivisionError):
        pass

    return out


def reporting_odds_ratio(table: ContingencyTable) -> dict:
    """ROR(报告比值比)及 95% 置信区间(对数正态近似)。

    公式:
        ROR    = (a / b) / (c / d) = a*d / (b*c)
        SE(lnROR) = sqrt(1/a + 1/b + 1/c + 1/d)
        95% CI = exp(ln ROR +/- 1.96 * SE)

    ROR 即 2x2 表的交叉乘积比,等价于该表优势比的估计量。
    出现零单元格时使用 Haldane-Anscombe(+0.5)校正后的表计算。
    """
    corrected = table.has_zero_cell
    t = table.corrected()

    out = {
        "ror": None,
        "ci_low": None,
        "ci_high": None,
        "se_ln_ror": None,
        "corrected": corrected,
    }

    # 空表不具备任何信息,返回 None 而不是校正后的 ROR=1
    if table.n <= 0 or t.b <= 0 or t.c <= 0:
        return out

    ror = (t.a * t.d) / (t.b * t.c)
    out["ror"] = ror

    try:
        se = math.sqrt(1.0 / t.a + 1.0 / t.b + 1.0 / t.c + 1.0 / t.d)
        out["se_ln_ror"] = se
        if ror > 0 and se > 0 and math.isfinite(se):
            ln = math.log(ror)
            out["ci_low"] = math.exp(ln - 1.96 * se)
            out["ci_high"] = math.exp(ln + 1.96 * se)
    except (ValueError, ZeroDivisionError):
        pass

    return out


def bcpnn_ic(table: ContingencyTable) -> dict:
    """BCPNN 的 Information Component(IC)及 IC025/IC975。

    采用 Bate 等(1998)提出的近似方差公式,即 WHO UMC 使用的标准做法:

        IC      = log2( a * N / ((a+b) * (a+c)) )         [即观测/期望 的 2 为底对数]

        E(IC)   = log2( (a + gamma11) * (N + alpha) * (N + beta)
                        / ((N + gamma) * (a + b + alpha1) * (a + c + beta1)) )

        V(IC)   = 1/(ln2)^2 * [ (N - a + gamma - gamma11) / ((a + gamma11) * (1 + N + gamma))
                                + (N - (a+b) + alpha - alpha1) / ((a + b + alpha1) * (1 + N + alpha))
                                + (N - (a+c) + beta  - beta1 ) / ((a + c + beta1 ) * (1 + N + beta )) ]

        IC025 = E(IC) - 1.96 * sqrt(V(IC))
        IC975 = E(IC) + 1.96 * sqrt(V(IC))

    先验参数取 UMC 常用的一组弱信息先验:
        gamma11 = 0.5, gamma = 1.0, alpha = 1.0, beta = 1.0, alpha1 = 0.5, beta1 = 0.5
    在 N 较大时,该式退化为 IC = log2(a / E),其中 E = (a+b)(a+c)/N。

    信号判据:IC025 > 0(即 IC 的 95% 下限为正)。
    """
    # BCPNN 先验超参数(UMC 惯例)
    gamma11, gamma, alpha, beta, alpha1, beta1 = 0.5, 1.0, 1.0, 1.0, 0.5, 0.5

    n = table.n
    out = {
        "ic": None,
        "ic_raw": None,
        "ic025": None,
        "ic975": None,
        "ic_var": None,
        "expected": None,
        "corrected": table.has_zero_cell,
    }

    if n <= 0 or table.a <= 0:
        # a=0 时观测/期望为 0,IC 无定义(理论上为 -inf)
        return out

    # 期望值用原始表计算,便于报告展示
    exp_count = table.n_drug * table.n_event / n
    out["expected"] = exp_count
    if exp_count > 0:
        out["ic_raw"] = math.log2(table.a / exp_count)

    a = table.a
    try:
        e_ic = math.log2(
            ((a + gamma11) * (n + alpha) * (n + beta))
            / ((n + gamma) * (a + table.b + alpha1) * (a + table.c + beta1))
        )
        var_ic = (1.0 / (math.log(2.0) ** 2)) * (
            (n - a + gamma - gamma11) / ((a + gamma11) * (1.0 + n + gamma))
            + (n - (a + table.b) + alpha - alpha1)
            / ((a + table.b + alpha1) * (1.0 + n + alpha))
            + (n - (a + table.c) + beta - beta1)
            / ((a + table.c + beta1) * (1.0 + n + beta))
        )
    except (ValueError, ZeroDivisionError):
        return out

    if var_ic < 0 or not math.isfinite(var_ic):
        return out

    sd = math.sqrt(var_ic)
    out["ic"] = e_ic
    out["ic_var"] = var_ic
    out["ic025"] = e_ic - 1.96 * sd
    out["ic975"] = e_ic + 1.96 * sd
    return out


# --------------------------------------------------------------------------
# 组合评估
# --------------------------------------------------------------------------
def contingency_table_from_counts(
    a: float,
    drug_total: float,
    event_total: float,
    grand_total: float,
) -> ContingencyTable:
    """由 a、药物总报告数、反应总报告数与全库总数构造 2x2 表。

    b = n1 - a  (目标药物的其他反应)
    c = m1 - a  (其他药物的目标反应)
    d = N - a - b - c
    任一结果出现负数说明输入口径不一致,统一以 max(0, .) 截断以避免静默产生错误的负数格。
    """
    b = max(0.0, drug_total - a)
    c = max(0.0, event_total - a)
    d = max(0.0, grand_total - a - b - c)
    return ContingencyTable(a=float(a), b=float(b), c=float(c), d=float(d))


def evaluate_pair(
    drug: str,
    reaction: str,
    table: ContingencyTable,
    require_all_criteria: bool = False,
) -> SignalResult:
    """对一个药物-反应 2x2 表计算全部统计量并套用信号判据。

    require_all_criteria=True 时,只有 PRR、ROR、IC 三条判据同时满足才算信号;
    默认 False,即任一判据满足即标记(但各判据的结果分别记录在 *_signal 字段)。
    """
    prr_res = proportional_reporting_ratio(table)
    ror_res = reporting_odds_ratio(table)
    ic_res = bcpnn_ic(table)

    res = SignalResult(
        drug=drug,
        reaction=reaction,
        a=table.a,
        b=table.b,
        c=table.c,
        d=table.d,
        n_total=table.n,
        prr=prr_res["prr"],
        prr_ci_low=prr_res["ci_low"],
        prr_ci_high=prr_res["ci_high"],
        chi2=prr_res["chi2"],
        chi2_p=prr_res["p_value"],
        ror=ror_res["ror"],
        ror_ci_low=ror_res["ci_low"],
        ror_ci_high=ror_res["ci_high"],
        ic=ic_res["ic"],
        ic025=ic_res["ic025"],
        ic975=ic_res["ic975"],
        ic_var=ic_res["ic_var"],
        expected=ic_res["expected"],
        corrected=bool(prr_res["corrected"]),
    )

    res.is_signal_prr = (
        res.prr is not None
        and res.chi2 is not None
        and res.prr >= EU_PRR_MIN
        and res.chi2 >= EU_CHI2_MIN
        and table.a >= EU_N_MIN
    )
    res.is_signal_ror = res.ror_ci_low is not None and res.ror_ci_low > 1.0
    res.is_signal_ic = res.ic025 is not None and res.ic025 > 0.0

    if require_all_criteria:
        keep = res.is_signal_prr and res.is_signal_ror and res.is_signal_ic
        res.is_signal_prr = res.is_signal_prr and keep
        res.is_signal_ror = res.is_signal_ror and keep
        res.is_signal_ic = res.is_signal_ic and keep

    return res


def sort_key_signal_strength(res: SignalResult) -> tuple:
    """按信号强度排序:优先 PRR,其次 IC025,最后报告数 a。降序。"""
    return (
        -(res.prr if res.prr is not None else 0.0),
        -(res.ic025 if res.ic025 is not None else -99.0),
        -res.a,
    )
