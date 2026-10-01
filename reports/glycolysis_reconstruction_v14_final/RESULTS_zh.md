# 无酶糖酵解代谢物网络重建：第一轮

共运行 30 个起点网络、153 次 seed→TS 尝试、348771 次势能/力评估，注册 30 条 AIMNet2-2025 TS 边；平均每条注册边 11626 次评估。

其中 20 条化学边、10 条构象边；7 条包含发起极小值，5 条严格命中符号提议端点。

## 已知片段审计

| 片段 | 证据类型 | 严格恢复 | 未完全收敛的精确端点候选 |
|---|---|---:|---:|
| G6P ↔ F6P | full_atom | 否 | 0 |
| FBP ↔ GAP + DHAP | full_atom | 是 | 1 |
| GAP ↔ DHAP | full_atom | 否 | 0 |
| 3PG ↔ 2PG | full_atom | 否 | 0 |
| 2PG ↔ PEP + H2O | full_atom | 是 | 6 |
| PEP + H2O ↔ pyruvate + H3PO4 | phosphate_surrogate | 否 | 0 |
| glucose + H3PO4 ↔ G6P + H2O | phosphate_surrogate | 否 | 4 |
| F6P + H3PO4 ↔ FBP + H2O | phosphate_surrogate | 否 | 2 |
| 1,3-BPG + H2O ↔ 3PG + H3PO4 | phosphate_surrogate | 否 | 1 |

严格恢复 2/9 个参考片段。

## 严格恢复的物理边

- **FBP ↔ GAP + DHAP**：双向候选能垒 2.754 / 2.264 eV；来源 `reports/glycolysis_targeted_refine_v14/fbp/cleavage_fbp_o0/arrows/network.json`。
- **2PG ↔ PEP + H2O**：双向候选能垒 1.900 / 2.043 eV；来源 `reports/glycolysis_targeted_refine_v14/dehydration/dehydration_2pg_o0/arrows/network.json`。

## 证据边界

- 所有物理边均为 AIMNet2-2025 鞍点与双侧下降，不是 DFT/IRC。
- 磷酸基团采用完全质子化中性微观状态；结果不代表生理 pH。
- H3PO4/H2O 只作为守恒原子的磷酸储库替代，不代表 ATP/ADP 酶促能垒。
- 未加入 NAD+/NADH，因此 GAPDH 氧化还原步骤不在本轮物理搜索中。