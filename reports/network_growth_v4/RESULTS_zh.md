# 扩展体系与多步网络增长 v4

参照 ReactionAtlas 的双来源提议、验证后扩展和新极小值队列设计；这是净中性小体系的受限实验，不是完整 formose、离子或溶液动力学复现。
参考：[ReactionAtlas 方法与补充实验](https://arxiv.org/html/2606.30778v1)。

## 实验设计

八个固定起点、几何/符号/混合三种策略，各一个随机种子。相同 16000 个几何评估预算，最多 24 次尝试；每个几何独立计费，批量调用不折算成一次评估。当前策略根据已观察网络分配动作，没有更新神经网络权重。

| 起点 | 策略 | 尝试 | 几何评估 | 后端调用 | 用时/s | 根连通物种 | 化学图对 | 最长已证路径/化学步 | 新物种续探 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| acetone | geometry | 9 | 10828 | 10538 | 107.6 | 1 | 0 | 0 | 0 |
| acetone | arrows | 18 | 16000 | 14376 | 172.4 | 2 | 1 | 1 | 1 |
| acetone | hybrid | 17 | 16000 | 14724 | 185.1 | 2 | 1 | 1 | 1 |
| propanal | geometry | 15 | 16000 | 15072 | 155.8 | 2 | 1 | 1 | 1 |
| propanal | arrows | 15 | 16000 | 14782 | 176.7 | 3 | 2 | 1 | 1 |
| propanal | hybrid | 15 | 16000 | 14956 | 171.4 | 3 | 3 | 2 | 2 |
| glycolaldehyde | geometry | 14 | 16000 | 15540 | 154.4 | 2 | 1 | 1 | 1 |
| glycolaldehyde | arrows | 18 | 16000 | 14221 | 174.3 | 4 | 4 | 2 | 1 |
| glycolaldehyde | hybrid | 17 | 16000 | 14620 | 170.5 | 3 | 3 | 1 | 2 |
| glyceraldehyde | geometry | 14 | 16000 | 15310 | 184.5 | 3 | 2 | 2 | 2 |
| glyceraldehyde | arrows | 13 | 16000 | 15655 | 206.4 | 2 | 1 | 1 | 1 |
| glyceraldehyde | hybrid | 13 | 16000 | 15655 | 208.0 | 2 | 1 | 1 | 1 |
| acetylacetone | geometry | 9 | 11495 | 11408 | 144.5 | 1 | 0 | 0 | 0 |
| acetylacetone | arrows | 14 | 16000 | 15913 | 223.9 | 1 | 0 | 0 | 0 |
| acetylacetone | hybrid | 14 | 16000 | 15913 | 222.2 | 1 | 0 | 0 | 0 |
| cyclobutanone | geometry | 14 | 16000 | 15370 | 187.6 | 2 | 1 | 1 | 1 |
| cyclobutanone | arrows | 17 | 16000 | 14268 | 190.8 | 4 | 3 | 2 | 3 |
| cyclobutanone | hybrid | 16 | 16000 | 14425 | 191.0 | 5 | 4 | 2 | 4 |
| acetamide | geometry | 9 | 10504 | 10400 | 115.9 | 1 | 0 | 0 | 0 |
| acetamide | arrows | 17 | 16000 | 14492 | 176.6 | 2 | 1 | 1 | 1 |
| acetamide | hybrid | 16 | 16000 | 14752 | 169.4 | 2 | 1 | 1 | 1 |
| nitroethane | geometry | 14 | 16000 | 15594 | 176.3 | 2 | 1 | 1 | 1 |
| nitroethane | arrows | 6 | 6399 | 5471 | 69.1 | 3 | 2 | 1 | 0 |
| nitroethane | hybrid | 15 | 16000 | 15014 | 169.9 | 4 | 3 | 2 | 3 |

物种包含立体信息；不同构象不重复算物种。两步增长要求实际极小值网络存在路径，且压缩连续构象后至少有三个不重复物种。不通过合并同 SMILES 的未连接构象制造通路，也不把 A→B→A 往返计为增长。

## 实际多步候选

- propanal / hybrid：CCC=O → [H]/[C-]=[O+]\CC → CC.[C-]#[O+]；实际节点 [0, 2, 3]；边 [1, 2]。
- glycolaldehyde / arrows：O=CCO → O/C=C\O → [H]/[O+]=C\C[O-]；实际节点 [0, 2, 9]；边 [1, 7]。
- glyceraldehyde / geometry：O=C[C@@H](O)CO → [H]/[O+]=C\[C@H]([O-])CO → [H]/[O+]=C\C([O-])CO；实际节点 [0, 1, 4]；边 [0, 2]。
- cyclobutanone / arrows：O=C1CCC1 → O=CC1CC1 → C1CC1.[C-]#[O+]；实际节点 [0, 2, 3]；边 [1, 2]。
- cyclobutanone / hybrid：O=C1CCC1 → O=CC1CC1 → C1CC1.[C-]#[O+]；实际节点 [0, 2, 3]；边 [1, 2]。
- nitroethane / hybrid：CC[N+](=O)[O-] → CCn1oo1 → CC/N=[O+]\[O-]；实际节点 [0, 2, 3]；边 [1, 2]。

## 独立 DFT

- propanal_edge1：physical_event_verified；端点对保留=False；['[H]/[C-]=[O+]/CC', 'CCC=O']。
- propanal_edge2：unresolved_ts；端点对保留=None；[]。
完整原链 DFT 认证：False；实际两步 DFT 链：False。
DFT 未收敛属于该计算尝试未认证，不是反应不存在的证据；替代端点与原端点分别保存。

## 主实验后的独立检查

下述检查不改写主实验尝试、预算或成功计数；当前代码已包含真实力停止条件的修正。

- 几何预算补查 acetamide：16000 个评估，根连通化学图对 0；原尝试前缀保持=True。
- 几何预算补查 acetone：16000 个评估，根连通化学图对 0；原尝试前缀保持=True。
- 几何预算补查 acetylacetone：16000 个评估，根连通化学图对 0；原尝试前缀保持=True。
- Hessian 三点检查通过=True，虚频数量保持；批量单项计时比约 5.91–6.94。并发负载会影响计时。
- Dimer 力外推 acetone：评估 799→685，状态 validated_descents→validated_descents，端点图对相同=True。
- Dimer 力外推 enol：评估 1079→995，状态 validated_descents→validated_descents，端点图对相同=True。
- Dimer 力外推 glycolaldehyde：评估 619→556，状态 validated_descents→validated_descents，端点图对相同=True。
- Dimer 力外推 cyclobutanone：评估 712→631，状态 validated_descents→validated_descents，端点图对相同=True。
- 停止条件检查 glycolaldehyde / 6：ts_force_unconverged→validated_descents，实际端点 ['O=CCO', 'O[C@H]1CO1']；同图端点不自动计为新化学连接。
- 停止条件检查 acetone / 1：ts_force_unconverged→validated_descents，实际端点 ['CC(C)=O', 'C=C(C)O']；同图端点不自动计为新化学连接。
- 停止条件检查 acetone / 10：ts_force_unconverged→validated_descents，实际端点 ['CCC=O', 'CCC=O']；同图端点不自动计为新化学连接。
- 修正后端到端短运行：8 次尝试、8000 次评估、4 条 MLIP 边。

## 效率与限制

总尝试 339；总几何评估 359226；后端调用 338469；整个并行活动耗时 1592.1 秒。单个运行用时包含算法、图处理和写盘，不含进程导入。
- Hessian 分批求力，保持原有限差分步长和认证阈值。计时不能外推成整体或 DFT 加速倍数。
- 使用预训练 AIMNet2-rxn member0 与 ASE Dimer；所有 TS/端点采用 0.005 eV/Å 力阈值和完整曲率检查。
- 网络相关优先级是启发式；混合组每四次中安排一次几何提议，并允许符号无覆盖节点的几何探索。
- 单随机种子不能支持统计显著性；新物种可能仍在预训练数据分布内。
- 没有把电子能垒当成自由能或实验反应速率；没有把搜索次数当成反应概率。
- 双向模型训练、独立箭头真值和反应物分子池的双分子增长仍属后续工作。

可视化：[总览](index.html)、[丙醛实际分子网络](search/propanal_s17/molecules/index.html)、[环丁酮实际分子网络](search/cyclobutanone_s17/molecules/index.html)。
