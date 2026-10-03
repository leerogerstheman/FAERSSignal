# FAERSSignal —— FAERS 药物警戒信号检测

基于美国 FDA **openFDA** 公开接口的真实不相称性分析(disproportionality analysis)工具。
对指定药物与不良事件组合,从零实现 **PRR / ROR / BCPNN-IC** 四项统计量及其置信区间,
输出 SQLite 数据库、CSV 明细与中文 Markdown 报告。

本项目不是演示用的假数据管线:仓库中的 `examples/` 目录保存了**一次真实运行**的
全部产物,README 中引用的每个数字都可以追溯到该次运行。

> **重要声明**:本工具是**信号筛查**工具,不产出临床结论,不证明因果关系,
> 不构成任何医疗建议。详见文末[局限性](#六局限性)。

---

## 一、真实运行结果(非示例编造)

以下数字来自 `examples/` 中保存的实际运行,分析对象为
**ASPIRIN、IBUPROFEN、METFORMIN** 三个药物 × 45 个 MedDRA 首选术语(PT),
共 135 个药物-反应组合。全库药物-反应报告对总数 **N = 39,497,182**。

### 汇总

| 指标 | 数值 |
|---|---|
| 药物-反应组合数 | 135 |
| 检出信号数(任一判据) | **125** |
| 三条判据(PRR+ROR+IC)同时满足 | **111** |
| 仅满足 PRR 判据 | 111 |
| 仅满足 ROR 判据 | 125 |
| 仅满足 IC025>0 判据 | 125 |

各药物分别的信号数:ASPIRIN 42/45,IBUPROFEN 41/45,METFORMIN 42/45。

### 最强信号

**METFORMIN —— 乳酸酸中毒(LACTIC ACIDOSIS)**,这是本次分析中不相称程度最高的组合:

- 报告数 **a = 20,130**(2x2 表为 a=20130, b=412767, c=11542, d=39052743)
- **PRR = 157.38**,95% CI 153.85–160.99
- **ROR = 165.01**,95% CI 161.24–168.86
- **IC = 5.86,IC025 = 5.83**,chi2 = 1,140,825.4
- 期望报告数 E = 347.13,即实际报告数为期望值的约 **58 倍**

### 其他代表性信号

| 药物 | 反应 (PT) | a | PRR | ROR (95% CI) | IC025 |
|---|---|---:|---:|---|---:|
| METFORMIN | LACTIC ACIDOSIS | 20,130 | 157.38 | 165.01 (161.24–168.86) | 5.83 |
| METFORMIN | HYPOGLYCAEMIA | 9,151 | 20.16 | 20.57 (20.11–21.05) | 4.03 |
| METFORMIN | BLOOD GLUCOSE INCREASED | 27,171 | 15.37 | 16.33 (16.12–16.55) | 3.71 |
| ASPIRIN | GASTROINTESTINAL HAEMORRHAGE | 14,998 | 10.34 | 10.57 (10.39–10.76) | 3.15 |
| IBUPROFEN | JOINT SWELLING | 6,491 | 10.37 | 10.65 (10.39–10.93) | 3.26 |
| IBUPROFEN | ACUTE KIDNEY INJURY | 6,889 | 8.63 | 8.87 (8.66–9.09) | 3.01 |
| ASPIRIN | MYOCARDIAL INFARCTION | 16,436 | 5.56 | 5.69 (5.60–5.78) | 2.35 |

这些结果与已知的临床认识方向一致(二甲双胍与乳酸酸中毒、NSAIDs 与消化道出血和肾损伤),
可以视为对统计管线正确性的一个旁证。但请注意:**这些数字反映的是报告不相称性,
不是风险大小,也不证明因果关系**——见[局限性](#六局限性)。

### 结果文件

| 文件 | 说明 |
|---|---|
| [`examples/signals.csv`](examples/signals.csv) | 135 行完整结果明细(含 2x2 表与全部统计量) |
| [`examples/report.md`](examples/report.md) | 自动生成的完整中文 Markdown 报告 |
| [`examples/faers_example.sqlite`](examples/faers_example.sqlite) | SQLite 数据库(analysis_run + signal_result 两张表) |

---

## 二、安装与使用

### 环境要求

- Python 3.9+(开发与验证环境为 Python 3.12)
- 仅依赖标准库即可运行分析;`pandas` 不是必需的
- 需要可访问 `https://api.fda.gov` 的网络

无需安装第三方依赖:

```powershell
# 项目根目录下直接运行
python faers_signal.py --drugs ASPIRIN,IBUPROFEN,METFORMIN --top-reactions 30 `
    --db examples/faers_example.sqlite `
    --csv examples/signals.csv `
    --md examples/report.md
```

### 复现 README 中的结果

上表数字对应的完整命令(在 `D:\FAERSSignal` 下执行):

```powershell
$env:PYTHONIOENCODING='utf-8'
python faers_signal.py --drugs ASPIRIN,IBUPROFEN,METFORMIN --top-reactions 30 `
    --db examples/faers_example.sqlite --csv examples/signals.csv --md examples/report.md
```

约需 **52 次 API 请求、95 秒**(默认 1 秒/请求的礼貌限流)。
由于 openFDA 数据滚动更新,后续复现的绝对数值会随报告累积而略有变化。

### 常用参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--drugs` | `ASPIRIN,IBUPROFEN,METFORMIN,ATORVASTATIN,SERTRALINE` | 逗号分隔的药物名(openFDA 规范名,大写) |
| `--top-reactions` | `30` | 每个药物取报告数最高的前 N 个反应作为候选 |
| `--db` | `data/faers.sqlite` | SQLite 输出路径 |
| `--csv` | `output/signals.csv` | CSV 输出路径 |
| `--md` | `output/report.md` | Markdown 报告路径 |
| `--sleep` | `1.0` | 两次请求之间的间隔秒数 |
| `--api-key` | 无 | 可选 openFDA API key,提高限流上限 |
| `--grand-total` | `0` | 直接指定全库总数 N,>0 时跳过该次查询 |
| `--top-n` | `40` | 报告中列出的信号条数上限 |

### 单元测试

```powershell
python tests/test_stats.py          # 独立运行,无需 pytest
python -m pytest tests -q           # 若已安装 pytest
```

当前状态:**27 个测试全部通过**。所有期望值均由公式手工推导,
并已用独立方法交叉验证(见[第五节](#五实现的验证))。

---

## 三、数据来源

数据来自 **openFDA `drug/event` 接口**,即美国 FDA 不良事件报告系统(FAERS)的公开数据。

- 接口文档:https://open.fda.gov/apis/drug/event/
- 计数通过 `count=` 聚合端点获得,例如:
  ```
  https://api.fda.gov/drug/event.json?search=patient.drug.medicinalproduct:"ASPIRIN"&count=patient.reaction.reactionmeddrapt.exact&limit=200
  ```

**数据口径说明(影响结果解读,必读)**

1. `a` 取自该药物的反应分布聚合结果中对应 PT 的计数。
2. `n1 = a+b`(该药物报告总数)与 `m1 = a+c`(该反应报告总数)取自
   `meta.results.total`,**不受聚合查询 `limit` 截断影响**。
3. `N`(全库药物-反应报告对总数)由按药物名聚合的前若干项求和得到。
   实测中 openFDA 对该聚合端点**拒绝 `limit=1000`(返回 HTTP 403)**,
   `limit<=800` 正常,故客户端按 800→500→200→100 的顺序降级尝试,
   本次运行实际使用 `limit=800`。由于未覆盖报告数极少的尾部药物名,
   **N 是下限估计**,这会使 IC 与卡方被轻微高估(PRR/ROR 的方向不受影响)。
4. 反应使用 MedDRA 首选术语(PT)层级的精确匹配,未合并 SMQ;
   药物名按 openFDA 规范化的 `medicinalproduct` 精确匹配,
   复方制剂、商品名与拼写差异可能造成漏计或串计。

---

## 四、方法学

### 2x2 列联表定义

对每一个"药物-反应"组合,按下表构造列联表:

```
                    目标反应 a        其他反应         合计
    目标药物           a               b            n1 = a+b
    其他药物           c               d            n0 = c+d
    合计            m1 = a+c        m0 = b+d         N = a+b+c+d
```

- `a`:该药物同时报告了该反应的报告数
- `b`:该药物报告了其他反应的报告数
- `c`:其他药物报告了该反应的报告数
- `d`:其余报告数

### 四项统计量

**1. PRR(比例报告比)**

```
PRR        = [ a/(a+b) ] / [ c/(c+d) ]
SE(ln PRR) = sqrt( 1/a - 1/(a+b) + 1/c - 1/(c+d) )
95% CI     = exp( ln PRR ± 1.96 · SE )
chi2       = N·(ad-bc)² / [ (a+b)(c+d)(a+c)(b+d) ]        (df = 1)
```

**2. ROR(报告比值比)**

```
ROR        = (a/b) / (c/d) = ad / (bc)
SE(ln ROR) = sqrt( 1/a + 1/b + 1/c + 1/d )
95% CI     = exp( ln ROR ± 1.96 · SE )
```

**3. BCPNN 的 IC(信息成分)**

```
IC      = log2( a / E ),      E = (a+b)(a+c) / N
E(IC)   = log2( (a+g11)(N+α)(N+β) / [ (N+γ)(a+b+α1)(a+c+β1) ] )
V(IC)   = 1/(ln2)² · [ (N-a+γ-g11)   / ((a+g11)(1+N+γ))
                       + (N-(a+b)+α-α1) / ((a+b+α1)(1+N+α))
                       + (N-(a+c)+β-β1) / ((a+c+β1)(1+N+β)) ]
IC025   = E(IC) - 1.96·sqrt(V(IC))
IC975   = E(IC) + 1.96·sqrt(V(IC))
```

先验超参数取 UMC 惯例的弱信息先验:`g11=0.5, γ=1.0, α=1.0, β=1.0, α1=0.5, β1=0.5`。
在 N 较大时该式退化为 `IC = log2(a/E)`。

**4. 卡方检验 p 值**

`chi2` 的 p 值由自实现的正则化下不完全伽马函数 `P(s,x)` 计算
(级数展开 + 连分式展开),未依赖 scipy。

### 信号判据

| 判据 | 阈值 | 出处 |
|---|---|---|
| PRR | PRR ≥ 2 且 chi2 ≥ 4 且 a ≥ 3 | 欧盟 EMA 常规信号阈值 |
| ROR | ROR 的 95% CI 下限 > 1 | 常规 ROR 判据 |
| IC | IC025 > 0 | WHO UMC BCPNN 判据 |

结果表中分别记录三条判据各自的命中情况,`PRR+ROR+IC` 表示三者同时满足。

### 连续性校正

当 2x2 表出现 0 单元格时,比值与对数会出现除零。此时对**四个单元格同时加 0.5**
(Haldane-Anscombe 校正),**仅在存在 0 单元格时启用**,相关行在报告中标注。
该校正使估计值向无效值(PRR=ROR=1)收缩,解读时需注意。

本次运行中有 4 行触发校正(均为 `a=0`,即该药物未报告该事件)。
此时 PRR/ROR 会得到极小的数值(如 0.001 量级),这是校正后的数学结果,
**不应解读为"该药物对此事件有保护作用"**——`a=0` 的真实含义是"未观察到报告"。

### 非特异性术语过滤

工具默认过滤 20 个反映报告流程或给药问题、而非医学不良事件的术语
(如 `DRUG INEFFECTIVE`、`OFF LABEL USE`、`PRODUCT QUALITY ISSUE`、
`MEDICATION ERROR` 等)。这类术语若纳入会产生大量无临床意义的"强信号"。

---

## 五、实现的验证

统计实现未依赖 scipy/statsmodels,因此对其正确性做了独立交叉验证:

| 验证项 | 独立参照方法 | 结果 |
|---|---|---|
| `chi2_sf(·, df=1)` | 卡方分布已知临界值 | 5 个临界点误差 < 1e-7 |
| `PRR` | Python `fractions` 精确有理数运算 | 3 组表格误差 < 1e-12 |
| `ROR` | 倒数对称性(交换行 ROR 取倒数) | 乘积 = 1,误差 < 1e-9 |
| `IC` | N 很大时应退化为 `log2(a/E)` | 差值 3.4e-4(量级符合预期) |
| `normal_ppf` | 二分法高精度参照 | 最大绝对误差 2.4e-7 |
| 全流程 | 重新调用 openFDA 原始计数手工复算 | 与流水线输出**完全一致** |

最后一项值得强调:对最强信号 METFORMIN→LACTIC ACIDOSIS,
我们从 openFDA 重新拉取原始计数(a=20130, n1=432897, m1=31672),
手工套用公式得到 PRR=157.383069、ROR=165.009626、chi2=1140825.4482,
与流水线输出的 2x2 表和全部统计量逐位吻合。

在开发过程中,上述验证也**发现并修复了两个真实缺陷**:
空表(四格全零)曾返回 `PRR=1`、`ROR=1` 这一具有误导性的"无关联"结论,
现改为返回 `None`;以及 `grand_total()` 的 `limit=1000` 触发 HTTP 403 的问题。

---

## 六、局限性

**这是解读任何结果的必要前提。**

1. **自发呈报数据没有分母**。FAERS 是自愿提交的报告集合,不存在暴露人群数,
   因此 PRR/ROR/IC **不是发病率、不是风险比、不是绝对风险**,
   也不能用于比较不同药物的安全性高低。
2. **不相称不等于因果**。结果仅表示"报告中出现了超出背景水平的聚集",
   不能证明药物导致了该事件。信号须经流行病学与临床评估才能形成安全性结论。
3. **报告偏倚**。报告行为受事件严重程度、上市时间、媒体关注、监管行动、
   诉讼与文献报道影响,新药与受关注药物更易被报告。
4. **重复报告**。同一病例可能由患者、医生、企业分别提交,FAERS 未完全去重,
   会同时抬高 `a` 与 `b`。
5. **适应症混杂**。药物所治疗的疾病本身可能导致该不良事件,产生虚假不相称。
6. **分组口径**。反应按 MedDRA PT 精确匹配,同一医学概念的不同 PT 会被拆散;
   药物名按规范化名称精确匹配,复方制剂与商品名可能漏计或串计。
7. **N 为下限估计**(见[第三节](#三数据来源)第 3 点),会使 IC 与卡方被轻微高估。
8. **零单元格校正的副作用**。`a=0` 的行经 +0.5 校正后会得到极小的 PRR/ROR,
   这是数学产物而非保护效应。
9. **多重比较**。对上百个组合同时检验会产生大量假阳性,
   本工具未做多重比较校正(药物警戒筛查阶段的惯例)。
10. **阈值是经验值**。PRR≥2 / chi2≥4 / IC025>0 等阈值是人为约定,
    不同机构的做法并不统一。
11. **数据时点**。openFDA 滚动更新,复现结果会随时间变化。

---

## 七、项目结构

```
D:\FAERSSignal\
├── faers_signal.py          入口脚本(插入项目根到 sys.path 后调用 src.cli.main)
├── src\
│   ├── __init__.py          包信息
│   ├── openfda.py           openFDA API 客户端(限流、退避重试、404 空结果处理)
│   ├── stats.py             PRR / ROR / BCPNN-IC 与自实现分布函数
│   ├── store.py             SQLite 持久化(analysis_run / signal_result)
│   ├── report.py            Markdown 报告生成
│   └── cli.py               命令行入口与流程编排
├── tests\
│   └── test_stats.py        27 个单元测试(手工推导期望值)
├── examples\                一次真实运行的产物
│   ├── signals.csv
│   ├── report.md
│   └── faers_example.sqlite
├── README.md
├── LICENSE                  MIT
└── .gitignore
```

---

## 八、参考文献

- Evans SJW, Waller PC, Davis S. Use of proportional reporting ratios (PRRs) for signal
  generation from spontaneous adverse drug reaction reports.
  *Pharmacoepidemiol Drug Saf*. 2001;10(6):483-486.
- Rothman KJ, Lanes S, Sacks ST. The reporting odds ratio and its advantages over the
  proportional reporting ratio. *Pharmacoepidemiol Drug Saf*. 2004;13(8):519-523.
- Bate A, Lindquist M, Edwards IR, et al. A Bayesian neural network method for adverse
  drug reaction signal generation. *Eur J Clin Pharmacol*. 1998;54(4):315-321.
- Norén GN, Hopstadius J, Bate A. Shrinkage observed-to-expected ratios for robust and
  transparent large-scale pattern discovery. *Stat Methods Med Res*. 2013;22(1):57-69.
- openFDA drug/event API:https://open.fda.gov/apis/drug/event/

---

## English Summary

**FAERSSignal** is a reproducible disproportionality-analysis toolkit for
pharmacovigilance, built on the public **openFDA `drug/event` API** (FAERS data).

It implements, from scratch and without scipy/statsmodels:

- **PRR** (Proportional Reporting Ratio) with chi-square, p-value and 95% CI
- **ROR** (Reporting Odds Ratio) with log-normal 95% CI
- **BCPNN IC** (Information Component) with the standard approximate variance,
  yielding IC025 / IC975
- Signal criteria: PRR≥2 & chi2≥4 & a≥3 (EU); ROR lower 95% CI > 1; IC025 > 0
- Haldane-Anscombe (+0.5) continuity correction, applied only when a zero cell is present

Outputs SQLite + CSV + a Chinese Markdown report (methods, contingency-table
definition, criteria, results sorted by signal strength, and limitations).

**Verified real run** (saved under `examples/`): 3 drugs × 45 MedDRA PTs = 135 pairs,
N = 39,497,182; **125 signals** detected, 111 satisfying all three criteria.
Strongest signal: **METFORMIN — LACTIC ACIDOSIS**, a = 20,130,
**PRR = 157.38** (95% CI 153.85–160.99), **ROR = 165.01** (95% CI 161.24–168.86),
**IC025 = 5.83**.

All statistics were cross-validated against independent references (exact rational
arithmetic, published chi-square critical values, and a full manual recomputation
from freshly fetched raw API counts). **27/27 unit tests pass.**

**Limitations**: FAERS is spontaneous reporting — there is no denominator, the
statistics are not incidence or risk, disproportionality is not causality, and the
data are subject to reporting bias, duplicate reports, indication confounding and
multiple comparisons. This tool is for signal screening only and is **not** medical
advice.

## License

MIT — see [LICENSE](LICENSE). Copyright (c) 2026 leerogerstheman.
