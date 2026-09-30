# 文献案例 benchmark v6

比较相同案例、相同模型与预算下的几何/净变键/完整箭头提议。FlowER 和 ReactionAtlas 原模型没有执行；本轮没有新增 DFT。
FlowER 小分子开环是论文 Fig.1b 的表示示例，不是单例 DFT 认证；ReactionAtlas 两例有论文中相应反应的讨论。
完整的 FlowER 十步示例含 Br/F/Cl，不能用当前势能；原测试集中符合严格大小/元素/组分筛选的一条离子对记录，在几何 Lewis 识别阶段受阻。失败与限制见 cases.json，未删原子简化。

| 案例 | 策略 | 命中 / 取向 | 尝试 | 几何评估 |
|---|---|---:|---:|---:|
| flower_epoxide_ammonia | 纯几何 | 0/2 | 14 | 18000 |
| flower_epoxide_ammonia | 净变键 | 0/2 | 16 | 18000 |
| flower_epoxide_ammonia | 完整箭头 | 0/2 | 17 | 18000 |
| atlas_glycolaldehyde_hydration | 纯几何 | 0/2 | 15 | 18000 |
| atlas_glycolaldehyde_hydration | 净变键 | 2/2 | 14 | 18000 |
| atlas_glycolaldehyde_hydration | 完整箭头 | 1/2 | 14 | 18000 |
| atlas_formaldehyde_dimer | 纯几何 | 0/2 | 15 | 18000 |
| atlas_formaldehyde_dimer | 净变键 | 0/2 | 17 | 18000 |
| atlas_formaldehyde_dimer | 完整箭头 | 0/2 | 16 | 18000 |

命中要求同一条通过 MLIP 检查的边同时匹配文献反应物与目标产物。仅出现目标分子、另一端不同，不算命中。
只有两个相遇取向，不能宣称统计泛化。完整箭头并未在本轮稳定优于净变键；环氧开环未通过体现当前 seed/搜索覆盖的不足。
不同轮次使用不同相遇构型和预算；v5 的甲醛二聚通过不能改写成 v6 命中。

总尝试 138；总几何评估 162000；主搜索墙钟 11.6 分钟。

统一入口：../index.html；文献与差距：../benchmark.html；实际分子视图：../events.html?stage=literature_benchmark_v6。
