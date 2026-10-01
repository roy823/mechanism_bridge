# AIMNet2-2025 正式算法反应路径扩展

`arrows` 表示完整箭头引导的三维 seed，不是只输出箭头。seed 使用净变键、活性原子、相遇取向、反应进度、进攻角以及电子源/受体耦合；随后由 AIMNet2-2025、Dimer 和双侧下降获得实际 B/C。

完成运行 13 组，237 次尝试，226512 个几何评估，登记 40 条 B–TS–C 边、24 个按体系区分的端点图对；其中 7 条边、5 个图对与选定文献案例相符。

## 查看保存的真实反应过程

- [羟基乙醛与甘油醛](intramolecular/molecules/index.html)
- [甲醛＋烯二醇与 Coley 环加成](targeted_bimolecular/molecules/index.html)
- [环氧乙烷＋氨](epoxide_ammonia/molecules/index.html)
- [甲醛＋羟基乙醛](bimolecular/molecules/index.html)

每个页面均使用保存的 TS 和双侧 BFGS 下降帧，没有插值。

## 文献相关命中

| 体系 | 实际端点 | 势垒/eV | 符号模板 | seed 预测端点命中 |
|---|---|---|---|---|
| formaldehyde_glycolaldehyde_o0 | `O=C[C@@H](O)CO` ↔ `C=O.O=CCO` | 4.264 / 3.658 | `grammar:carbonyl_C_H_coupled_addition` | 是 |
| flower_epoxide_ammonia_o1 | `NCCO` ↔ `C1CO1.N` | 3.198 / 2.196 | `grammar:epoxide_amine_opening_stepwise` | 否 |
| glycolaldehyde | `C=O.C=O` ↔ `O=CCO` | 1.899 / 4.482 | `SynEPD:4:reverse` | 否 |
| glycolaldehyde | `C=O.C=O` ↔ `O=CCO` | 1.899 / 4.482 | `SynEPD:4:reverse` | 否 |
| coley_3217_o1 | `CN1CCC=N1` ↔ `C#[N+][N-]C.C=C` | 2.854 / 0.385 | `grammar:dipolar_3_plus_2_cycloaddition` | 是 |
| formaldehyde_enediol_o1 | `O=C[C@@H](O)CO` ↔ `C=O.O/C=C/O` | 1.339 / 0.407 | `grammar:neutral_enol_aldol_addition` | 是 |
| formaldehyde_enediol_o1 | `O=C[C@@H](O)CO` ↔ `C=O.O/C=C/O` | 1.338 / 0.407 | `grammar:neutral_enol_aldol_addition` | 是 |

## 证据边界

- 新边是 AIMNet2-2025 势能面上的一阶鞍点和双侧极小值，不是 DFT/IRC。
- `source_connected` 仅说明实际端点是否包含发起节点 A；不作为 B–TS–C 事件的拒绝条件。
- 相同图对允许多条不同 TS；高能分解和替代通道继续保留。

## 案例来源

- 甲醛、羟基乙醛、烯二醇和甘油醛来自 ReactionAtlas 讨论的 formose 化学空间。
- 环氧乙烷＋氨来自 FlowER 的表示示例；该单例不是 FlowER 作者逐例 DFT 认证。
- `coley_3217` 来自 Coley 组公开的偶极 [3+2] 环加成数据，并有参考 TS 可用于后续独立复核。