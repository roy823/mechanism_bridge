# 双分子探索边界验证 v5

## 方法和例子来源

ReactionAtlas 提供甲醛、水、小糖及其相遇复合物的研究背景；采用其随机刚体取向、接近与验证后扩展的思路。本实验使用 AIMNet2-rxn 与 Dimer，没有运行 MoreRed/MD-ET，也没有复现溶液 formose 动力学。[原文](https://arxiv.org/html/2606.30778v1)。

Coley 组实际数据选取 rxn_id 3217、3216：按中性闭壳层 CHNO 的显式原子数、反应 ID 排序，取前两个不同反应物组合。参考产物仅在搜索完成后比较。[数据与代码](https://github.com/coleygroup/dipolar_cycloaddition_dataset)。
之前讨论的 Jin/Coley 强化学习论文案例为自由基氧化与自由基串联环化，超出本轮模型适用范围。[MIT 论文](https://hdl.handle.net/1721.1/151666)。

六组反应物、两种相对取向、三个策略；每次运行最多 12 次尝试和 12000 个几何评估。每个尝试最多 2000 次评估，物理力收敛阈值 0.005 eV/Å，TS 恰一个低于 -30 cm⁻¹ 的虚频，两端零显著虚频。初始极小值失败也保留。

符号动作包括公开 SynEPD 局部迁移与明确标注的人为电子动作规则；后者不是模型预测或独立金标准。符号组包含额外的反应位点取向采样，因此本轮比较整个符号引导方案，不单独归因于完整箭头。

## 实际结果

主实验共 250 次尝试、302414 个几何评估，墙钟 20.6 分钟。

| 起点 | 策略 | 状态 | 尝试 | 跨分子重原子事件 | 原始反应物匹配 | 严格根连通 | 评估 |
|---|---|---|---:|---:|---:|---:|---:|
| formaldehyde_water_o0 | geometry | completed | 9 | 0 | 0 | 0 | 11644 |
| formaldehyde_water_o0 | arrows | completed | 6 | 1 | 1 | 0 | 6210 |
| formaldehyde_water_o0 | hybrid | completed | 11 | 1 | 1 | 0 | 12000 |
| formaldehyde_water_o1 | geometry | completed | 9 | 0 | 0 | 0 | 11359 |
| formaldehyde_water_o1 | arrows | completed | 6 | 1 | 1 | 0 | 6557 |
| formaldehyde_water_o1 | hybrid | completed | 11 | 1 | 1 | 0 | 12000 |
| formaldehyde_dimer_o0 | geometry | completed | 9 | 0 | 0 | 0 | 11100 |
| formaldehyde_dimer_o0 | arrows | completed | 11 | 2 | 2 | 0 | 12000 |
| formaldehyde_dimer_o0 | hybrid | completed | 11 | 1 | 1 | 0 | 12000 |
| formaldehyde_dimer_o1 | geometry | initial_minimum_unresolved | 0 | 0 | 0 | 0 | 426 |
| formaldehyde_dimer_o1 | arrows | initial_minimum_unresolved | 0 | 0 | 0 | 0 | 426 |
| formaldehyde_dimer_o1 | hybrid | initial_minimum_unresolved | 0 | 0 | 0 | 0 | 426 |
| formaldehyde_glycolaldehyde_o0 | geometry | completed | 10 | 0 | 0 | 0 | 12000 |
| formaldehyde_glycolaldehyde_o0 | arrows | completed | 10 | 1 | 1 | 0 | 12000 |
| formaldehyde_glycolaldehyde_o0 | hybrid | completed | 10 | 1 | 1 | 0 | 12000 |
| formaldehyde_glycolaldehyde_o1 | geometry | completed | 9 | 0 | 0 | 0 | 11099 |
| formaldehyde_glycolaldehyde_o1 | arrows | completed | 10 | 0 | 0 | 0 | 12000 |
| formaldehyde_glycolaldehyde_o1 | hybrid | completed | 10 | 0 | 0 | 0 | 12000 |
| formaldehyde_enediol_o0 | geometry | completed | 9 | 0 | 0 | 0 | 11605 |
| formaldehyde_enediol_o0 | arrows | completed | 10 | 0 | 0 | 0 | 12000 |
| formaldehyde_enediol_o0 | hybrid | completed | 10 | 0 | 0 | 0 | 12000 |
| formaldehyde_enediol_o1 | geometry | completed | 9 | 0 | 0 | 0 | 11661 |
| formaldehyde_enediol_o1 | arrows | completed | 10 | 2 | 2 | 0 | 12000 |
| formaldehyde_enediol_o1 | hybrid | completed | 10 | 2 | 2 | 0 | 12000 |
| coley_3217_o0 | geometry | initial_minimum_unresolved | 0 | 0 | 0 | 0 | 335 |
| coley_3217_o0 | arrows | initial_minimum_unresolved | 0 | 0 | 0 | 0 | 335 |
| coley_3217_o0 | hybrid | initial_minimum_unresolved | 0 | 0 | 0 | 0 | 335 |
| coley_3217_o1 | geometry | completed | 9 | 0 | 0 | 0 | 11629 |
| coley_3217_o1 | arrows | completed | 6 | 1 | 1 | 0 | 7491 |
| coley_3217_o1 | hybrid | completed | 10 | 1 | 1 | 0 | 12000 |
| coley_3216_o0 | geometry | initial_minimum_unresolved | 0 | 0 | 0 | 0 | 341 |
| coley_3216_o0 | arrows | initial_minimum_unresolved | 0 | 0 | 0 | 0 | 341 |
| coley_3216_o0 | hybrid | initial_minimum_unresolved | 0 | 0 | 0 | 0 | 341 |
| coley_3216_o1 | geometry | completed | 9 | 0 | 0 | 0 | 11128 |
| coley_3216_o1 | arrows | completed | 6 | 0 | 0 | 0 | 7625 |
| coley_3216_o1 | hybrid | completed | 10 | 0 | 0 | 0 | 12000 |

“原始反应物匹配”允许已枚举的共振等价，但不把同 SMILES 的不同相遇构型直接合并为同一极小值。“严格根连通”要求真实极小值图上的通路；脱离根的发现单独保留。重原子成键、氢转移、旁观分子存在时的单分子重排分开统计。
事件数可包含同一化学图对的多个 TS，不能当作不同反应家族数。含孤立原子/离子、超过两个片段或高形式电荷的候选在结构化结果中另有标记；通过 MLIP 数值检查不代表它们具有可靠化学意义。

## 跨策略汇总

| 策略 | 尝试 | 几何评估 | 跨分子重原子事件 | 成功运行 / 12 | 严格根连通事件 |
|---|---:|---:|---:|---:|---:|
| geometry | 82 | 104327 | 0 | 0 | 0 |
| arrows | 75 | 88985 | 8 | 6 | 0 |
| hybrid | 93 | 109102 | 7 | 6 | 0 |

跨所有取向和策略去重后有 6 个重原子成键图对；不等于相互独立的反应家族。

## 实际反应图对

- `C#[N+][N-]C.C=C` ↔ `C=CN(C)N=C`
- `C=O.C=O` ↔ `O=CCO`
- `C=O.O` ↔ `OCO`
- `C=O.O/C=C/O` ↔ `O/C=C/OCO`
- `C=O.O/C=C/O` ↔ `O=C[C@H](O)CO`
- `C=O.O=CCO` ↔ `COC(=O)CO`

## 独立 DFT 检查

| 检查 | 状态 | 物理事件通过 | 原图对保留 | 端点 | 气相电子能垒/eV |
|---|---|---|---|---|---|
| formaldehyde_dimer | unresolved_ts | False | None |  |  |
| formaldehyde_dimer_refined | physical_event_verified | True | True | O=CCO / C=O.C=O | 4.5931, 3.8329 |
| formaldehyde_water | physical_event_verified | True | True | OCO / C=O.O | 2.2210, 1.5317 |

水合为预先选择的主验证；二聚为事后补充，严格优化另存。能垒按上表端点顺序，参照实际相遇极小值，不是无限分离的单体能量之和。
已验证事件坐标、Hessian 频率、端点与提议来源导出到 `DFT_events.jsonl`；失败记录保留于 `DFT_feedback.jsonl`。箭头仍是提议假设，没有独立金标准标签。

## Coley 参考 TS 诊断（不计入主实验）

| 条件 | 记录 | 结果 | 参考图对保留 | 评估 |
|---|---|---|---|---:|
| coley_reference_ts_diagnostic | 3217 | unresolved_minimum | True | 1060 |
| coley_reference_ts_diagnostic | 3216 | unresolved_minimum | True | 886 |
| coley_reference_ts_strict | 3217 | validated_descents | True | 1251 |
| coley_reference_ts_strict | 3216 | unresolved_minimum | True | 1485 |

两组都使用作者参考 TS 和虚频方向；strict 将真实力阈值从 0.005 收紧到 0.001 eV/Å。3217 的通过证明该局部 MLIP 通路可以被保留，不能算作我们的反应物起始搜索命中。3216 仍有端点优化问题。

## 相遇盆地诊断

首次 BFGSLineSearch/0.0005 的五组检查均未同时收敛，作为失败诊断保留。随后 BFGS/0.001 的结果：

| 体系 | 两端均为合格极小值 | RMSD/Å | 能量差/eV | 合并条件满足 |
|---|---|---:|---:|---|
| formaldehyde_water | True | 1.063 | 0.07543 | False |
| formaldehyde_dimer | True | 0.711 | 0.05515 | False |
| formaldehyde_glycolaldehyde | True | 1.767 | 0.02021 | False |
| formaldehyde_enediol | True | 0.787 | 0.02142 | False |
| coley_3217 | True | 0.984 | 0.03870 | False |

这些诊断没有修改主实验网络。下一步需要实际探索相遇复合物之间的构象转换或分离/再相遇过程，不能用较宽 RMSD 阈值代替物理连接。

## 结论边界

气相 MLIP 电子能垒不能直接与 Coley 数据表中的水相活化自由能比较。下降轨迹不是 IRC；只有明确列出的独立 DFT 检查可升级证据。没有神经模型更新、浓度动力学或产率预测，也没有声称穷举全部符号空间。

每次运行仍固定原子库存；本轮验证两个分子之间的探索，没有动态加入第三个分子，也不把断裂片段自动拼接成已验证网络边。
