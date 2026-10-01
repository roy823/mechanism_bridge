# Formose 最后缺口：中性 tetrose ↔ 2GO

本轮只使用 AIMNet2‑2025。没有运行 AIMNet2‑rxn 或 ensemble 交叉筛选。

在 11 个构象起点上完成 48 次定向任务和 50875 次势能/力评估，登记 5 条边；其中 1 条严格连接中性 tetrose 与 2GO。

## 严格命中

`O=CCO.O=CCO` ↔ `O=C[C@@H](O)[C@H](O)CO`，双向势垒 2.482 / 2.867 eV。符号提议端点命中：是。

命中只出现在 4 个 threose 构象中的一个；4 个 erythrose 构象和 3 个 2GO 逆向构象均没有严格命中。

## 循环结论

- **组成级：** 结合 v12 的前五段，现在 6/6 分段均有物理边，得到完整的 GO 自催化增殖候选。
- **立体分辨级：** 尚未闭合。C4 生长段到 D‑erythrose，本轮闭环边来自 D‑threose，二者之间没有物理互变边。

## 查看

- [全部 C2/C3/C4 分子和 92 条边的 TransitionNet](global_transitionnet/molecules/index.html)
- [命中闭环 TS 与真实双侧下降](runs/threose_uff_c3/molecules/index.html)

## 证据边界

- 当前是 AIMNet2‑2025 气相 MLIP 一阶鞍点与双侧下降，不是 DFT/IRC。
- 跨库存总网不比较节点绝对能量；各事件页面保留各自双向相对势垒。