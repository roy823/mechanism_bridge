# FlowER 论文反应类型的 GPA 式物理搜索

本轮选择 FlowER 论文明确展示的 Diels–Alder、Prilezhaev 环氧化和转酰胺化。使用的是三类反应的最小代表体系，不是论文受许可限制的逐条 NameRxn 测试分子。

流程为：完整电子箭头 → 净变键与活性原子 → 反应物相遇取向和三维 seed → AIMNet2-2025/Dimer → 一阶鞍点检查 → 双侧下降 → 登记实际 B–TS–C。搜索没有读取参考 TS 或产物几何。

共完成 6 个取向运行、84 次尝试和 83319 次势能/力评估；登记 13 条 TS 边，其中 5 条命中各自符号提议，4 条属于三类案例的目标事件。

## 结果

| 反应类型 | 代表反应 | 尝试/评估 | TS 边 | 结论 |
|---|---|---:|---:|---|
| Diels-Alder cycloaddition | `C=CC=C.C=C` | 24 / 25933 | 7 | 发现 1 条反应物↔环己烯目标边，并保留开链与带电旁路。 |
| Prilezhaev epoxidation | `C=C.OOC=O` | 12 / 13008 | 3 | 两个取向发现 2 条烯烃＋过酸↔环氧化物＋甲酸目标边。 |
| Transamidation | `CC(=O)N.CN` | 48 / 44378 | 3 | 找到胺加成的中性四面体事件；未形成完整取代产物与三步网络。 |

## 目标物理事件

| 类型 | 实际 B ↔ C | 双向势垒/eV | 符号模板 |
|---|---|---:|---|
| Diels-Alder cycloaddition | `C=C.C=CC=C` ↔ `C1=CCCCC1` | 1.467 / 0.993 | `grammar:diels_alder_4_plus_2_cycloaddition` |
| Prilezhaev epoxidation | `C=C.O=COO` ↔ `C1CO1.O=CO` | 1.738 / 3.982 | `grammar:prilezhaev_epoxidation` |
| Prilezhaev epoxidation | `C=C.O=COO` ↔ `C1CO1.O=CO` | 1.738 / 4.029 | `grammar:prilezhaev_epoxidation` |
| Transamidation | `CC(N)=O.CN` ↔ `CN[C@@](C)(N)O` | 1.915 / 1.257 | `SynEPD:907:forward` |

## 查看真实构型和 TransitionNet

- [物种聚合总网](aggregate/molecules/index.html)
- [Diels–Alder 全部 TS 与下降轨迹](diels_alder/molecules/index.html)
- [Prilezhaev 全部 TS 与下降轨迹](prilezhaev/molecules/index.html)
- [转酰胺化全部 TS 与下降轨迹](transamidation/molecules/index.html)

## 当前证据能说明什么

- 三种 FlowER 符号反应类型都至少有一个符号引导 seed 落到了化学相关的一阶鞍点事件。
- Diels–Alder 与 Prilezhaev 的完整目标端点对已出现；说明多箭头协同动作可以转成有效三维 seed。
- 转酰胺化只得到加成/互变相关局部事件，尚未自动形成完整的加成、质子转移和离去网络。它直接暴露了下一步应做的“新节点同步扩展”。
- 所有能垒均为 AIMNet2-2025 气相电子能结果；双侧 BFGS 下降不是严格 IRC，也没有完成 DFT 复核。