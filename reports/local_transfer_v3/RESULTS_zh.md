# 局部电子动作迁移：v3 工程验证

固定起点为丙醛、丙酮，各一个随机种子。所有对照在同一起点使用相同势、位移幅度与预算。
本轮改变了局部提议、seed 角度/异步进度和初始方向，不能与历史 v2 直接比较归因。

## 符号覆盖

固定面板 22 个分子中 16 个有合法提议。该计数只表示图级回放合法，不证明存在相应 TS。
记录审计：{'source_records': 1926, 'outside_pilot_domain': 1850, 'replay_valid_records': 48, 'replay_or_parse_rejected': 27, 'unsupported_local_directions': 11, 'no_distinct_chemical_graph': 1}。来源文件和全部动作见 coverage.json。

## 实际搜索

| 起点 | 策略 | 尝试 | 势评估 | 初始化评估 | 根连通化学图对 | 从新物种继续 |
|---|---|---:|---:|---:|---:|---|
| propanal | geometry | 6 | 4746 | 145 | 0 | False |
| propanal | center_random | 6 | 5279 | 145 | 0 | False |
| propanal | bond_edits | 6 | 5299 | 145 | 0 | False |
| propanal | arrows | 6 | 5248 | 145 | 0 | False |
| acetone | geometry | 6 | 5002 | 192 | 1 | False |
| acetone | center_random | 6 | 5079 | 192 | 0 | False |
| acetone | bond_edits | 6 | 3844 | 192 | 1 | False |
| acetone | arrows | 6 | 4232 | 192 | 1 | False |

图对按运行统计，反向事件和构象不算新的独立化学体系。候选数不是实际机理准确率。

## 已观察连接

- propanal / geometry：无根连通化学连接；结果分类 {'not_index_one': 4, 'calculation_failed': 1, 'same_basin_return': 1}。
- propanal / center_random：无根连通化学连接；结果分类 {'not_index_one': 6}。
- propanal / bond_edits：无根连通化学连接；结果分类 {'not_index_one': 4, 'new_connection': 1, 'unresolved_minimum': 1}。
- propanal / arrows：无根连通化学连接；结果分类 {'not_index_one': 4, 'unresolved_minimum': 2}。
- acetone / geometry：CC(C)=O ↔ CCC=O；结果分类 {'unresolved_minimum': 1, 'not_index_one': 3, 'calculation_failed': 1, 'new_connection': 1}。
- acetone / center_random：无根连通化学连接；结果分类 {'not_index_one': 6}。
- acetone / bond_edits：C=C(C)O ↔ CC(C)=O；结果分类 {'not_index_one': 2, 'unresolved_minimum': 1, 'calculation_failed': 1, 'new_connection': 2}。
- acetone / arrows：C=C(C)O ↔ CC(C)=O；结果分类 {'not_index_one': 3, 'calculation_failed': 1, 'new_connection': 2}。

## DFT

- acetone_arrows_edge0_DFT：physical_event_verified；物理连接认证=True；原 MLIP 端点对保留=True；实际端点=['CC(C)=O', 'C=C(C)O']。
  方法 wb97x/6-31g(d)，gas_phase；两侧电子能垒 [3.2030376179418454, 2.432796020252681] eV。
DFT 与 AIMNet2-rxn 训练参考层级不同；能垒差不能全部归因于模型误差。DFT 复核没有产生独立电子箭头真值。

## 独立续探诊断

直接使用新发现的烯醇端点，继续 3 次箭头搜索，消耗 1495 次势评估。
结果：{'not_index_one': 3}；实际连接：[]。
该探针在主试验结束后安排，不加入主对照，不作为独立反应体系或长期网络覆盖证据。

### 续探负模诊断

前三个拒绝候选重新计算 Hessian，共 186 次势评估，不改变原始判定。
前两例的主负模约 −2110 cm⁻¹，与初始 Cartesian 方向的绝对重合约 0.83；额外负模约 −66～−79 cm⁻¹，主要位移在另一侧甲基氢。
这是反应方向已接近、但横向柔性运动未稳定的线索；仅频率和方向重合不能独立完成模式归属。下一步可检验保持反应方向的横向曲率优化，不能直接忽略额外负模。
含该诊断的 MLIP 总评估数为 41938。

### 收紧阈值的单独验证

选择第一个被拒绝的续探候选，将力阈值收紧为 0.005 eV/Å，重新 Dimer 与双侧下降，消耗 644 次势评估。
结果 validated_descents；实际端点 ['CC(C)=O', 'C=C(C)O']；新端点与续探源极小值的几何/能量匹配见 continuation_refinement/connection_audit.json。
该结果支持这一例的失败与收敛精度有关。它是原可逆图对的后续恢复，不是新的独立化学反应，也没有独立 DFT 验证这一个后处理候选。
主对照 48 次和初次续探 3 次的原始判定均不修改。

## 起点诊断与证据边界

- 首轮 search 中两个起点均未通过极小值检查，原始输出保留。
- search_refined 收紧起点力阈值；丙醛通过，丙酮仍有两个负模。
- search_curvature 对丙酮沿负模微扰后重新最小化，额外计算计入初始化预算。
- 两体系的最终运行路径及源码快照分别保存在 summary.json 指向的目录中。
- 只有一个随机种子、两个相关体系；不作统计显著性或未见反应家族泛化声明。
- 局部模式仅支持当前非芳香活动中心和中性闭壳层 CHNO 域；有限提议上限会截断动作空间。
- 严格初始电子源检查排除了部分原库的中继式表达；不等于判定这些反应不可能。
- 完整箭头仍是提议条件，未成为独立审核的电子机理真值；未更新模型权重。

主对照共 48 次尝试、38729 次势评估；连同起点诊断共 39613 次势评估。

真实坐标可视化：[分子结构、TS 和下降轨迹](search_curvature/acetone_s17/molecules/index.html)。
