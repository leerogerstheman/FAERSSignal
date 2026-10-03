# FAERSSignal v1.0.0 — FAERS 药物警戒信号检测

基于 openFDA 全库的**不相称性分析**工具，从零实现 PRR / ROR / BCPNN 三种方法，
输出可复现的信号检测报告。**纯 Python 标准库，不依赖 scipy / statsmodels。**

这不是一个 ICSR 录入界面，而是一次真实的安全性分析 —— 从数据里算出信号。

## 实测结果（真实数据）

对 3 个药物 × 45 个首选术语（PT）= **135 个药物-反应组合**，
全库药物-反应报告对总数 **N = 39,497,182**：

| 指标 | 数值 |
| --- | --- |
| 检出信号数（任一判据） | **125 / 135** |
| 三条判据同时满足 | **111** |
| METFORMIN | 42 / 45 |
| IBUPROFEN | 41 / 45 |
| ASPIRIN | 42 / 45 |

**最强信号：METFORMIN → LACTIC ACIDOSIS**

| 统计量 | 数值 |
| --- | --- |
| a | 20,130 |
| **PRR** | **157.38**（95% CI 153.85–160.99） |
| **ROR** | **165.01**（95% CI 161.24–168.86） |
| **IC025** | **5.83** |
| χ² | 1,140,825.4 |
| 期望值 E | 347.13（实际约为期望的 58 倍） |

其他代表性信号：METFORMIN/HYPOGLYCAEMIA PRR=20.16；
ASPIRIN/GASTROINTESTINAL HAEMORRHAGE PRR=10.34；
IBUPROFEN/JOINT SWELLING PRR=10.37。

> 方向与已知临床认识一致（二甲双胍—乳酸酸中毒、NSAIDs—消化道出血/肾损伤），
> 构成对分析管线正确性的旁证。重跑两次结果完全一致。

## 实现的统计方法

全部从零实现，README 中给出每个公式：

- **PRR** — 比例报告比，含 χ² 与 95% CI
- **ROR** — 报告比值比，95% CI 用对数正态近似
- **BCPNN** — 贝叶斯置信传播神经网络，实现 IC 与 IC025 及近似方差
- **判据** — EU 标准（PRR≥2 且 χ²≥4 且 N≥3）、ROR 下限>1、IC025>0
- **零值处理** — Haldane-Anscombe 连续性校正，空表返回 `None`

## 验证强度

- **27 项单元测试**，用手工推导的期望值断言（非代码自生成快照）
- `chi2_sf` 对 5 个已知临界值误差 < 1e-7
- `PRR` 对 `fractions` 精确有理数误差 < 1e-12
- `IC` 在大 N 下退化为 `log2(a/E)`，符合理论（差 3.4e-4）
- 重新拉取原始 API 计数**手工复算**最强信号，与管线输出逐位吻合

## 修复的真实缺陷

1. **`grand_total()` 用 `limit=1000` 时 openFDA 返回 HTTP 403** — 硬阻塞，不修则
   整个分析无法运行。实测 `limit≤800` 正常，改为 800→500→200→100 降级重试。
2. **空表经 +0.5 校正后返回 PRR=1/ROR=1** — 这是误导性的"无关联"结论，改为返回 `None`。
3. **`EU_SIGNAL_CRITERIA` 在 `__all__` 中但未定义** — `import *` 会 AttributeError。

## 快速开始

```bash
git clone https://github.com/leerogerstheman/FAERSSignal
cd FAERSSignal

python faers_signal.py --drugs ASPIRIN,IBUPROFEN,METFORMIN --top-reactions 30 \
  --db examples/faers_example.sqlite --csv examples/signals.csv --md examples/report.md
```

约 52 次请求 / 95 秒（1 秒礼貌限流）。

```bash
python tests/test_stats.py     # 27/27 passed
```

## 局限（重要）

- FAERS 是**自发报告**系统，**没有分母**，不能计算发生率
- 报告数受报告偏倚、媒体关注、诉讼影响；存在重复报告与信息不全
- **不相称性不等于因果关系**，信号仅提示需进一步评估
- N 为下限估计（仅取报告数最高的前 800 个药物名求和），会使 IC/χ² 轻微高估
- 未做多重比较校正，因此高命中率属该方法的预期行为（筛查工具而非因果推断）
- 4 行因 a=0 触发连续性校正得到 0.001 量级的 PRR/ROR，这是校正的数学产物，
  **不能**解读为保护作用
- 数据滚动更新，后续复现绝对值会略有变化

## 许可

MIT
