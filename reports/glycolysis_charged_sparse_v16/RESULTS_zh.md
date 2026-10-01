# 带电稀疏糖酵解骨架搜索

- 初始 roots：12（6 个化学起点，各 2 个构象/取向）
- 完成网络：18（含 6 个 recovery/continuation 网络）
- AIMNet2-2025 势能/力评估：380610
- 原始物理极小值/TS 边：89/47
- 聚合物种节点/边：19/27
- 实际端点匹配符号提议的事件：3/47

## 目标通路
- 通过：G6P <-> F6P skeleton connection；direct=True，shortest_path=1 edges
- 通过：GAP <-> DHAP skeleton connection；direct=True，shortest_path=1 edges
- 通过：F6P + H2PO4- <-> FBP + H2O；direct=False，shortest_path=2 edges
- 未通过：stable G6P-enediol-F6P endpoint sequence；逐步边=[]
- 未通过：stable GAP-enediol-DHAP endpoint sequence；逐步边=[False, True]

## 解释边界

- 每个双分子 seed 只含一个化学计量 Pi 或一个水分子，没有溶剂浴。
- TS 要求全体系收敛且仅有一个显著虚频；端点允许按碳骨架力阈值验收。
- 这是 AIMNet2-2025 气相 MLIP 双侧下降证据，不是酶催化机理或 DFT IRC。
- 聚合 HTML 为每条聚合边保留最多 61 个采样构型帧，并单独保存 TS XYZ。
- 旧版 F6P 输入的磷酸位置错误，旧 G6P–F6P 阴性结果不再作为证据。
