# 稀疏端点糖酵解网络重建

- 初始化学 roots：9；构象/取向 starts：18
- 完成局部网络：32
- AIMNet2-2025 势能/力评估：635834
- 原始 TS 边：117
- 聚合固定库存网络：7；物种/边：93/115
- 隐藏目标覆盖：9/9

## 九段审计

- 通过：G6P ↔ F6P；direct=True；path=1 edges
- 通过：FBP ↔ GAP + DHAP；direct=False；path=2 edges
- 通过：GAP ↔ DHAP；direct=True；path=1 edges
- 通过：3PG ↔ 2PG；direct=False；path=2 edges
- 通过：2PG ↔ PEP + H2O；direct=True；path=1 edges
- 通过：PEP + H2O ↔ pyruvate + Pi；direct=False；path=4 edges
- 通过：glucose + Pi ↔ G6P + H2O；direct=True；path=1 edges
- 通过：FBP + H2O ↔ F6P + Pi；direct=False；path=2 edges
- 通过：1,3-BPG + H2O ↔ 3PG + Pi；direct=False；path=2 edges

## 证据边界

- 搜索只读取单侧 root、符号动作和实际发现节点；隐藏目标仅用于最终审计。
- 物理证据为 AIMNet2-2025 气相一阶鞍点与双侧下降，不是酶催化自由能或 DFT IRC。
- ATP/ADP、NAD+/NADH、Mg2+ 和蛋白环境尚未加入；磷酸化步骤是骨架替代反应。