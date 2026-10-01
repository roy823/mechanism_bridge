# 无酶糖酵解代谢物网络重建：第一轮

共运行 52 个起点网络、282 次 seed→TS 尝试、672005 次势能/力评估，注册 86 条 AIMNet2-2025 TS 边；平均每条注册边 7814 次评估。

其中 67 条化学边、19 条构象边；12 条包含发起极小值，23 条严格命中符号提议端点，4 条仅按含碳骨架验收。

## 已知片段审计

| 片段 | 证据类型 | 直接边 | 多步可达 | 未完全收敛的精确端点候选 |
|---|---|---:|---:|---:|
| G6P ↔ F6P | full_atom | 否 | 否 | 0 |
| FBP ↔ GAP + DHAP | full_atom | 是 | 是 | 1 |
| GAP ↔ DHAP | full_atom | 否 | 是 | 0 |
| 3PG ↔ 2PG | full_atom | 是 | 是 | 0 |
| 2PG ↔ PEP + H2O | full_atom | 是 | 是 | 6 |
| PEP + H2O ↔ pyruvate + H3PO4 | phosphate_surrogate | 是 | 是 | 1 |
| glucose + H3PO4 ↔ G6P + H2O | phosphate_surrogate | 是 | 是 | 4 |
| F6P + H3PO4 ↔ FBP + H2O | phosphate_surrogate | 否 | 否 | 3 |
| 1,3-BPG + H2O ↔ 3PG + H3PO4 | phosphate_surrogate | 是 | 是 | 1 |

严格恢复 7/9 个参考片段。

## 严格恢复的物理边

- **FBP ↔ GAP + DHAP**：双向候选能垒 2.754 / 2.264 eV；端点验收 `legacy_full_system`；来源 `reports/glycolysis_targeted_refine_v14/fbp/cleavage_fbp_o0/arrows/network.json`。
- **3PG ↔ 2PG**：双向候选能垒 2.149 / 2.354 eV；端点验收 `validated_descents`；来源 `reports/glycolysis_relaxed_recovery_v15/runs/pg_shift_2pg_o0/pg_shift_2pg_o0/arrows/network.json`。
- **3PG ↔ 2PG**：双向候选能垒 2.516 / 2.543 eV；端点验收 `validated_descents`；来源 `reports/glycolysis_relaxed_recovery_v15/runs/pg_shift_3pg_o0/pg_shift_3pg_o0/arrows/network.json`。
- **3PG ↔ 2PG**：双向候选能垒 2.377 / 2.373 eV；端点验收 `validated_descents`；来源 `reports/glycolysis_relaxed_recovery_v15/runs/pg_shift_3pg_o0/pg_shift_3pg_o0/arrows/network.json`。
- **2PG ↔ PEP + H2O**：双向候选能垒 1.900 / 2.043 eV；端点验收 `legacy_full_system`；来源 `reports/glycolysis_targeted_refine_v14/dehydration/dehydration_2pg_o0/arrows/network.json`。
- **PEP + H2O ↔ pyruvate + H3PO4**：双向候选能垒 2.499 / 2.145 eV；端点验收 `validated_descents`；来源 `reports/glycolysis_relaxed_recovery_v15/runs/pep_hydrolysis_pep_water_o1/pep_hydrolysis_pep_water_o1/arrows/network.json`。
- **glucose + H3PO4 ↔ G6P + H2O**：双向候选能垒 2.482 / 2.498 eV；端点验收 `validated_descents`；来源 `reports/glycolysis_relaxed_recovery_v15/runs/glucose_phosphorylation_g6p_water_o0/glucose_phosphorylation_g6p_water_o0/arrows/network.json`。
- **glucose + H3PO4 ↔ G6P + H2O**：双向候选能垒 2.359 / 2.491 eV；端点验收 `validated_descents`；来源 `reports/glycolysis_relaxed_recovery_v15/runs/glucose_phosphorylation_g6p_water_o0/glucose_phosphorylation_g6p_water_o0/arrows/network.json`。
- **glucose + H3PO4 ↔ G6P + H2O**：双向候选能垒 2.322 / 2.339 eV；端点验收 `validated_descents`；来源 `reports/glycolysis_relaxed_recovery_v15/runs/glucose_phosphorylation_g6p_water_o1/glucose_phosphorylation_g6p_water_o1/arrows/network.json`。
- **glucose + H3PO4 ↔ G6P + H2O**：双向候选能垒 2.321 / 2.455 eV；端点验收 `validated_descents`；来源 `reports/glycolysis_relaxed_recovery_v15/runs/glucose_phosphorylation_g6p_water_o1/glucose_phosphorylation_g6p_water_o1/arrows/network.json`。
- **1,3-BPG + H2O ↔ 3PG + H3PO4**：双向候选能垒 1.619 / 1.817 eV；端点验收 `validated_descents`；来源 `reports/glycolysis_relaxed_recovery_v15_resume/bpg_water_o1/bpg_hydrolysis_bpg_water_o1/arrows/network.json`。

## 证据边界

- 所有物理边均为 AIMNet2-2025 鞍点与双侧下降，不是 DFT/IRC。
- 磷酸基团采用完全质子化中性微观状态；结果不代表生理 pH。
- H3PO4/H2O 只作为守恒原子的磷酸储库替代，不代表 ATP/ADP 酶促能垒。
- 未加入 NAD+/NADH，因此 GAPDH 氧化还原步骤不在本轮物理搜索中。