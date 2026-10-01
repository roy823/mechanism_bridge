# 并行共享 C3H6O3 formose / 三碳糖 TransitionNet

完整细胞糖酵解包含十个酶促步骤以及 ATP/ADP、NAD⁺/NADH、无机磷酸、Mg²⁺和带电磷酸糖，不适合作为当前中性 CHNO 气相原型的直接测试。本轮选择 ReactionAtlas 同样关注的 formose 化学空间，在固定 C3H6O3 原子库存上检验较长连续网络。

使用 4 个 CPU worker（每个 5 线程），完成 120 个唯一 reservation 和 115009 次势能/力评估。得到 44 个物理极小值、39 个构象簇、26 个物种和 30 条 TS 边。

## 网络拓扑

物种图有 3 个连通分量；起点所在分量覆盖 22/26 个物种，最远最短路径为 10 步，并出现 2 个简单回路。

30 条边中 29 条由完整箭头 seed 找到、1 条由几何 seed 找到；15 条边的实际端点命中发起符号提议。这个计数没有同预算纯几何对照，不能单独解释为因果增益。

## 一条最深物种路径

| 步 | A | B | 双向势垒/eV | seed |
|---:|---|---|---:|---|
| 1 | `C=O.O/C=C/O` | `O=C[C@@H](O)CO` | 1.338 / 0.407 | `grammar:neutral_enol_aldol_addition` |
| 2 | `O=C[C@@H](O)CO` | `O.O=C/C=C/O` | 2.855 / 3.178 | `SynEPD:812:reverse` |
| 3 | `O.O=C/C=C/O` | `O/C=C/C(O)O` | 1.885 / 1.586 | `SynEPD:270:forward` |
| 4 | `O/C=C/C(O)O` | `O=CCC(O)O` | 2.256 / 2.708 | `SynEPD:4:reverse` |
| 5 | `O=CCC(O)O` | `C=CO.O=CO` | 1.126 / 1.016 | `geometry` |
| 6 | `C=CO.O=CO` | `CC=O.O=CO` | 0.897 / 0.487 | `grammar:carbonyl_addition_with_proton_transfer` |
| 7 | `CC=O.O=CO` | `C[C@H](O)C(=O)O` | 3.082 / 2.477 | `grammar:carbonyl_C_H_coupled_addition` |
| 8 | `C[C@H](O)C(=O)O` | `CC(O)=C(O)O` | 1.837 / 2.988 | `SynEPD:17:forward` |
| 9 | `CC(O)=C(O)O` | `C[C@@H](O)C(=O)O` | 2.449 / 3.509 | `SynEPD:17:reverse` |
| 10 | `C[C@@H](O)C(=O)O` | `CC=O.[H]/[O+]=[C-]/O` | 2.187 / 0.101 | `SynEPD:4:reverse` |

第一步为甲醛＋烯二醇与甘油醛之间的羟醛成键事件。后续包含水合/脱水、羰基–烯醇互变、裂解和重组。

## 查看

- [真实 TS、双侧下降、物种合并 TransitionNet](molecules/index.html)

## 解释边界

- 这是无磷酸、无酶的 formose/三碳糖气相子网，不是糖酵解十步生物通路。
- “10 步路径”位于物种合并图；同一物种的不同构象簇之间不一定已有显式构象转换 TS。
- 所有新边均为 AIMNet2-2025 证据，没有进行 DFT/IRC 复核。