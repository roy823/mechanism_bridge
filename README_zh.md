# Mechanism Bridge

**统一结果入口：[研究总览 HTML](reports/index.html)**。查看 [全部物理事件与真实三维轨迹](reports/events.html)、[benchmark 方向、文献差距与同题案例](reports/benchmark.html)。重建页面运行 `python scripts/exploration/build_report_site.py`；说明见 [统一报告与 benchmark](docs/BENCHMARK_AND_REPORTS_zh.md)。

整站离线文件与本轮计算已打包：[证据包与恢复说明](reports/repository_snapshots/report_portal_v6_README.md)。

双分子探索 v5：[方法与边界](docs/BIMOLECULAR_V5_zh.md)、[实际结果](reports/bimolecular_v5/RESULTS_zh.md)、[分子与三维路径](reports/bimolecular_v5/index.html)。扩展至 ReactionAtlas 相关 CHO 反应物组合和 Coley 数据中的双分子环加成；区分真实跨分子成键、旁观物重排与未连通的相遇构型。

下一阶段设计：[符号动作、GPA 式事件网络与 MCTS](docs/MCTS_REACTION_SEARCH_DESIGN_zh.md)。该文档是待实现方案，包含搜索定义、能量评分、动力学产品判定和验证计划。

文献与实现衔接：[ReactionAtlas 数据核查、符号到 seed 与 GPA 式闭环](docs/REACTIONATLAS_SYMBOLIC_GPA_zh.md)。包含公开数据版本差异、事件证据分层和下一阶段验证设计。

符号与几何引导的分子反应网络探索：从当前极小值提出种子，用机器学习势寻找鞍点，以实际双侧下降端点登记连接，再继续探索新节点。对代表事件补独立 DFT/IRC 验证。

## 当前状态

- **已经运行：** RGD1 起始构型、SynEPD 箭头回放、预训练 AIMNet2-rxn、Dimer、双侧下降、实际端点网络与三维可视化。
- **验证 v2：** 位移幅度统一、反应中心控制、随机重复、扩展起点及独立 DFT/IRC。
- **局部迁移 v3：** 局部价态/电荷 SMARTS、初始电子源占据检查、异步距离与角度 seed、seed 处的内坐标初始方向；实际输出见 [v3 报告](reports/local_transfer_v3/RESULTS_zh.md)。
- **尚未实现：** 训练好的逆向箭头模型、独立审核的大规模箭头配对集、FlowER/RitS/React-OT 接入。

| 结果 | 入口 |
|---|---|
| 网络探索 v1：三起点、99 次尝试 | [报告](reports/network_exploration_v1/RESULTS_zh.md)、[分子视图](reports/network_exploration_v1/molecules/index.html) |
| 验证 v2 | [固定方案](docs/VALIDATION_V2_zh.md)、[实际结果](reports/validation_v2/RESULTS_zh.md) |
| 局部迁移 v3：丙醛与丙酮，48 次主对照及单独续探诊断 | [实施说明](docs/LOCAL_TRANSFER_V3_zh.md)、[结果与 DFT 证据](reports/local_transfer_v3/RESULTS_zh.md)、[真实分子视图](reports/local_transfer_v3/search_curvature/acetone_s17/molecules/index.html) |
| 网络增长 v4：8 体系、24 组运行及效率检查 | [方案](docs/NETWORK_GROWTH_V4_zh.md)、[实际结果](reports/network_growth_v4/RESULTS_zh.md) |
| 历史 RF、DFT 和 IBO 可行性试验 | [历史结果](FEASIBILITY_RESULTS_zh.md) |

v1 中符号组优于当时的随机几何基线，但位移幅度未统一，且完整箭头尚无稳定超越净变键的证据。v2 专门检验这些限制。乙醛与乙烯醇属于同一反应家族，不算独立样本。

## 仓库结构

```text
src/mechbridge/       可复用的数据、提议、搜索、验证和绘图模块
scripts/data/         下载、导入和起点准备
scripts/exploration/  网络探索、配置化试验、汇总和可视化
scripts/qc/           DFT/IRC 与轨道分析
scripts/baselines/    历史图学习和反馈探针
scripts/diagnostics/  数据及运行环境诊断
configs/             固定的试验设计
tests/               回归检查
assets/              离线三维查看器和许可证
data/                原始及处理数据
models/              固定版本权重
reports/             实际输出、来源、源码快照
```

详见 [架构](docs/ARCHITECTURE.md)。整理前源码保存在 `reports/repository_snapshots/pre_validation_v2.zip`；原始数据与历史结果未覆盖。旧入口已移动，未保留兼容脚本。新实验各自保存源码 ZIP、输入和模型哈希。

## 环境

Python 3.11+。Windows 探索环境与 Linux/WSL 量化环境分别固定依赖。

```powershell
python -m venv .venv-explore
.venv-explore\Scripts\python.exe -m pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu
.venv-explore\Scripts\python.exe -m pip install -r requirements-exploration-lock.txt
.venv-explore\Scripts\python.exe -m pip install -e . --no-deps
```

Linux/WSL 使用 `.venv-qc` 和 `requirements-qc.txt`。当前探索为 CPU 上的 AIMNet2-rxn 单个成员；DFT 标准为气相闭壳层 ωB97X/6-31G(d)。两者能量分开报告。

## 常用命令

已有原始数据时准备起点（下载来源见 [数据说明](DATA_DOWNLOADS_zh.md)）：

```powershell
.venv-explore\Scripts\python.exe scripts/data/fetch_exploration_resources.py
.venv-explore\Scripts\python.exe scripts/data/prepare_feasibility.py --max-atoms 8 --limit 12
.venv-explore\Scripts\python.exe scripts/data/prepare_network_pilot.py
.venv-explore\Scripts\python.exe scripts/data/prepare_validation_starts.py
```

运行固定对照方案，输出目录必须是新的：

```powershell
.venv-explore\Scripts\python.exe scripts/exploration/run_validation_campaign.py --config configs/validation_v2.json --outdir reports/validation_v2/search
```

独立运行与绘图：

```powershell
.venv-explore\Scripts\python.exe scripts/exploration/run_network_exploration.py --outdir reports/new_run --starts data/processed/validation_starts.jsonl --start-ids MR_8342_1 --strategies geometry center_random bond_edits arrows --attempts 6 --seeds-per-node 3 --attempt-evaluations 1200
.venv-explore\Scripts\python.exe scripts/exploration/summarize_network_exploration.py reports/new_run
.venv-explore\Scripts\python.exe scripts/exploration/render_molecular_exploration.py reports/new_run
```

独立 DFT/IRC（WSL 内）：

```bash
.venv-qc/bin/python scripts/qc/verify_network_connection.py --run reports/network_exploration_v1 --start MR_8342_1 --strategy arrows --edge 0 --outdir reports/new_dft_verification --selection "One predefined representative acetaldehyde connection"
```

## 证据等级

1. MLIP 一阶鞍点：局部驻点与曲率检查。
2. MLIP 双侧下降连接：两侧负模位移后到达极小值，不是严格 IRC。
3. DFT/IRC 连接：指定 DFT 标准上的精修、路径和端点证据。
4. 电子结构支持的箭头假设：需另做轨道或响应分析及稳定性检查。
5. 独立审核的箭头配对：当前还没有大规模真值集。

由 A 提出的种子若连接 B/C，只登记 B/C。失败、未收敛、符号库不覆盖与意外有效连接各自记录。电子守恒不能代替机理准确率。

首轮限 CHNO、净中性闭壳层；不覆盖溶液、自由基和电子态交叉。AIMNet2-rxn 预训练包含 RGD1，当前试验不证明未见化学泛化。当前符号库使用局部电荷/价态模式，排除芳香活动中心与初始电子源不满足的表达；覆盖不到的新物种会停止符号分支。

更多资料：[探索方法](docs/GPA_REACTION_EXPLORATION_zh.md)、[历史复现](docs/FEASIBILITY_RUN_zh.md)、[文献重评](docs/LITERATURE_REASSESSMENT_2026-09-29_zh.md)、[TS 到机理方案](docs/TS_TO_MECHANISM_WORKFLOW_zh.md)。
