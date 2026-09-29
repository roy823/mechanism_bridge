# 可行性实验复现

本次实验以原理验证为目标，不依赖 GPA 的具体实现。图层基线与量化证据分开记录，不将净键变化标成弯箭头真值，不将自动轨道解释当成人工标签。

## 数据与环境

核心数据来源及校验见 `../DATA_DOWNLOADS_zh.md` 与 `../reports/local_download_inventory.json`。原始 RGD1 有 176,992 条组记录；Transition1x v4 实际读取到 10,073 个反应组，官方 train/val/test 分别为 9,561/225/287。`data` 为全集，不能再与三个 split 拼接。

额外获取 [Tsuneda & Taketsugu, Communications Chemistry 2025](https://doi.org/10.1038/s42004-025-01556-5) 的官方补充 PDF 与原始作图数据，保存在 `data/raw/roef_2025/`，下载来源与 SHA256 见其 `receipt.json`。它们是文献对照资料，未作为箭头训练标签；完整逐帧数据在论文中注明需联系作者。本次未发送申请。

另检索了 [PMechDB 下载页](https://deeprxn.ics.uci.edu/pmechdb/download) 和 [质子转移数据 v1](https://doi.org/10.6084/m9.figshare.30875087.v1)。PMechDB 的表单尚未提交；质子转移包是约 943 MB 的符号资源，缺少本实验需要的对应 TS/IRC，因此仅保存元数据，没有追加整包下载。

计算环境位于项目内 `.venv-qc`，运行于已有 WSL Ubuntu。Windows 全局 Python 未替换。依赖版本见 `requirements-qc.txt`，下载的 Linux wheel 在 `.cache/qc-wheels`，哈希在 `reports/qc_wheel_sha256.json`。WSL 网络不可用，本次通过 Windows 下载 wheel，再在 WSL 离线安装。

从项目目录运行以下命令（可进入 WSL 后执行，也可使用 `wsl.exe -d Ubuntu --cd '/mnt/d/学习相关文档/LEARN/MolGPA/mechanism_bridge' --` 前缀）：

普通联网 Linux 可先执行 `python3.12 -m venv .venv-qc`，再用 `.venv-qc/bin/pip install -r requirements-qc.txt`。本机已有环境无需重装；断网 WSL 使用 `.venv-qc/bin/python tools/pip.pyz install --no-index --find-links .cache/qc-wheels -r requirements-qc.txt`。

```bash
.venv-qc/bin/python -m pytest -q --basetemp work/pytest-new-run
.venv-qc/bin/python scripts/data/prepare_feasibility.py --max-atoms 8 --limit 12
.venv-qc/bin/python scripts/baselines/benchmark_graph_bridge.py
.venv-qc/bin/python scripts/qc/run_feasibility.py --limit 3 --threads 2 --outdir reports/new-quantum-run
```

最后一条计算会占用 CPU，具体时间取决于 SCF、TS 和 IRC 收敛。已有量化结果不会覆盖；新实验使用新的输出目录。后处理脚本当前读取正式运行目录 `reports/feasibility/`。

## 数据筛选

`prepare_feasibility.py` 从真实 HDF5 选择最多 8 原子、单片段、净中性且 Lewis 图闭壳层的事件。检查两端三维结构感知得到的图与源 SMILES 一致，排除恒等端点、分数键级等情况。映射号定义为三维原子索引加一，采用源坐标身份；后续真正计算的路径保持原子索引。

输出 `feasibility_pool.jsonl` 用于图层基线；`feasibility_events.jsonl` 按原子数、事件 ID 确定小批量子任务，选择不参考计算结果。单重态是明确的试验假设，未声称完成电子态的全局或自旋稳定性认证。相同组成的源记录在同一折，逆反应、构象不会跨折。

## 图层双向模型

`benchmark_graph_bridge.py` 使用随机森林建立两个代理任务：

- 正向：反应物图/几何特征，加上或不加净键变化，预测源 TS 原子对距离。
- 反向：反应物特征，加上或不加 TS 原子对距离，预测带符号的键级净变化。

五折验证按分子组成隔离。`heldout_predictions.npz` 是留出结果；`graph_proxy_models.joblib` 是随后使用全部图层数据训练的模型，不用于生成上述留出指标。寻鞍点测试使用留出预测，而不是全数据模型。

正向比较中加入净键变化等于增加产物图信息，不能证明完整电子流优于已有端点图。反向目标是净键变化，不能充当弯箭头准确率。226 条源记录中有 192 个不同端点图对，不能把 5,501 个原子对当成独立反应样本；误差差值采用分子组成组 bootstrap。

## 物理事件与电子结构

`verification.py` 在同一 ωB97X/6-31G(d) 势上优化 TS、计算解析 Hessian、投影分子刚体模态，并用 [Sella IRC](https://github.com/zadorlab/sella/wiki/IRC) 沿正反两向积分。IRC 内迭代失败即记录未解决，不降级为最小化后仍称 IRC。两端再优化、检查最小值并与预期图匹配。

`electronic.py` 计算沿路径的 IAO/IBO，并按交叉 AO 重叠跟踪轨道。离散电子结构样本包括 IRC 与已优化端点，明确保存采样索引。`event_graph.py` 将电子对分配给端点 Lewis 位点容量，产生守恒的源—汇箭头假设；这里的电子对守恒由构造保证，不是独立验证指标。

低轨道重叠需要审阅。可运行：

```bash
.venv-qc/bin/python scripts/qc/audit_orbitals.py --event MR_453777_0
```

该脚本比较单轨道和整个占据子空间的重叠，并检查最后一步的近等分数匹配。所保存的候选并不穷尽所有等价表示，也不宣称其中某个就是唯一机理。

还可运行 `scripts/qc/check_representation_stability.py --event MR_453777_0`，将整条路径做相同刚体变换后重新计算 IBO，以测试原始箭头集合的稳定性。若集合改变，需要区分有限 DFT 网格误差、轨道局域化和跟踪的影响，不能仅凭守恒认为标签已经唯一确定。

`contract_atom_relays` 实现一个待评估的受限代数约定：仅收缩一入一出且电子对净变化为零的原子位点，保留中继记录，不合并不同电子对源—汇对应或消除闭环。本次旋转前后的 4 支/2 支箭头在该约定下得到相同流对应；这不构成一般化学等价性的证明。

## 留出模型初猜的真实计算

```bash
.venv-qc/bin/python scripts/baselines/benchmark_seed_recovery.py --event MR_453777_0 --budget 60
```

模型由留出预测的 TS 距离矩阵经三维嵌入得到初猜；对照为对齐后的端点插值。两者都不使用参考 TS 作初猜，各有 60 次 DFT 梯度预算。解析 Hessian 的附加验证成本单独存在，不能只凭梯度次数给出总体算力加速比。事后按驻点、阶数、参考结构 RMSD 和能量差判断是否找回已验证鞍点；这项恢复指标不等于重新执行了每个种子的完整 IRC。

初猜对照和刚体变换试验拒绝覆盖已有结果；再次运行时使用 `--outdir` 指定新目录。正式报告默认汇总原运行目录。

## 输出与结论

```bash
.venv-qc/bin/python scripts/qc/summarize_feasibility.py
```

自动生成 `FEASIBILITY_RESULTS_zh.md`、`reports/feasibility_results.json`、图层对照图和逐事件路径证据图。报告区分：图层互补信息、物理连接验证、自动符号假设、独立箭头真值、多轮反馈学习。后两项不能由本次小样本直接证明。

`data/processed/verified_event_pilot.jsonl` 保存以物理验证事件为单位的试验数据。不同路径实例用独立 ID 保存，相同化学端点组通过 `chemical_event_group` 联系；箭头尚未独立审核，因此 `verified_pair` 保持 false。

`data/processed/physics_feedback_records.jsonl` 保存物理结果对原查询的支持或不相容证据。原查询与新事件不相容时，新物理事件仍保留；该记录只是查询—事件配对的负证据，不是“原反应不可能”的标签，也不宣称发现了文献中未知的新反应。

模型找到了不同一阶鞍点时，可以继续检查其实际连接，而不是直接标成失败反应：

```bash
.venv-qc/bin/python scripts/qc/verify_alternative.py --parent-event MR_375774_0
```

新证据更新探针：

```bash
.venv-qc/bin/python scripts/baselines/probe_physics_feedback.py
```

该脚本在至少三个独立组成组具备新量化证据时，执行一次按组成留一验证的几何残差更新。源 B3LYP-D3/TZVP 几何模型与新 ωB97X/6-31G(d) 标签分开标注；它是跨层级校准的小样本探针，不是完整主动学习实验。意外连接不会被当作原目标的几何标签加入。

性能候选试验 `check_density_reuse.py` 比较冷初猜与复用上一构型密度。三点检查中力差异约 1.08e-4 eV/Å，超过预设 1e-5 阈值，因此正式实验保持 `reuse_density=False`。结果保留在 `reports/density_reuse_check.json`，没有为加速而放宽物理判据。
