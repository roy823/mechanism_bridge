# Mechanism Bridge v0.1

本机下载请以 [DATA_DOWNLOADS_zh.md](DATA_DOWNLOADS_zh.md) 为准；当前状态在 `reports/local_download_inventory.json`。代码可行性、文献及模型清单见 [本机评估](docs/FEASIBILITY_AND_MODELS_zh.md)。下文量化运行结果属于原始打包记录。

本机新增的真实数据试验见 [运行结果](FEASIBILITY_RESULTS_zh.md) 与 [复现说明](docs/FEASIBILITY_RUN_zh.md)：包括图层双向基线、DFT/IRC 连接验证、IBO 箭头假设、留出初猜寻鞍点对照，以及意外事件的后续验证。它们与下文历史运行记录分别保存。

后续主线见 [TS 到机理的模型与算法方案](docs/TS_TO_MECHANISM_WORKFLOW_zh.md)：随机森林仅作基线，重点转向现有反应势、TS 电子响应和事件级符号反推。

把符号机理与三维反应事件放到同一套可追溯的数据结构中，为后续“候选生成 → 鞍点验证 → 电子结构解释 → 配对学习”做准备。

**这是一套已经运行过真实数据导入和量化计算的研究原型；目前没有训练好的箭头预测模型，也没有完成高可信箭头—TS 配对集。** 当前不使用反应模板生成机理，不自动补试剂，不把键级差称作弯箭头真值。

## 已实现内容

| 模块 | 文件 | 实际行为 |
|---|---|---|
| 下载与校验 | `scripts/download_core.py` | 固定五个核心文件；aria2 并行续传；大小、MD5、SHA256 核对；独立本机完成清单 |
| 数据导入 | `src/mechbridge/adapters.py` | mech-USPTO CSV、FlowER 文本、RGD1 HDF5/CSV、Transition1x HDF5 |
| 守恒检查 | `chemistry.py` | 组成、电荷、映射编号唯一性、映射原子身份；保留显式氢和立体信息 |
| 候选配对 | `pairing.py` | 完整反应体系的端点图匹配，允许正逆向；输出候选，不升级为物理配对 |
| 数据划分 | `pairing.py` | 相同端点对、逆反应、构象、同一 sequence 通过连通分量放在同侧 |
| 三维搜索/检查 | `physics.py` | ASE Dimer；有限差分 Hessian；质量加权投影平动/转动；负模态双侧 BFGS |
| 计算后端 | `backends.py` | PySCF 闭壳层气相能量/力；自定义 ASE MLIP 工厂接口 |
| 电子结构 | `electronic.py` | IAO 基底中的四次方 Pipek–Mezey/IBO；轨道布居；相邻构型重叠与 Hungarian 对应 |
| 训练输入 | `training_data.py` | 流式读取 Transition1x 能量/力，不把优化采样误作 IRC |
| 计算任务 | `scripts/prepare_jobs.py` | 输出可审核的 XYZ、计算参数及命令；不隐含指定电荷和自旋 |
| 动力学基础 | `kinetics.py` | 显式活化自由能的单分子 TST 与速率矩阵；拒绝电子能垒输入 |

## 安装

推荐 Python 3.12、Linux 或 WSL2。基础数据处理不需要 GPU。

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[qc,test]'
pytest -q
```

`requirements-tested.txt` 记录原始打包环境版本，本机环境见 `reports/local_code_audit.json`。Windows 应在 WSL2 运行 PySCF。初始代码包包含每库 1,000 条处理示例、历史计算输出和脚本；原始数据需要单独下载，模型权重尚未安装。

## 数据下载与放置

预期大小、哈希和来源见 `reports/data_inventory.json` 与 `docs/DATASETS.md`；本机实际完成状态见 `reports/local_download_inventory.json`。原始文件放在：

```text
data/raw/mech_uspto/mech-USPTO-31k.csv
data/raw/flower/data.zip
data/raw/rgd1_zenodo/RGD1_allrxns.h5
data/raw/rgd1_zenodo/RGD1CHNO_AMsmiles.csv
data/raw/transition1x/Transition1x.h5
```

从作者公开接口直接获取完整文件，共 **8.30 GB**，无需临时附件的 gzip 分卷。FlowER ZIP 包含旧/新两个数据变体，可直接流式读取，不必全部展开（展开约 5.97 GB）。

```bash
python scripts/download_core.py --connections 4 --workers 3
```

上述命令使用项目内 Windows aria2 工具；Linux 使用 `--aria2 /path/to/aria2c`。中断后运行同一条命令续传，保留 `.part` 和 `.part.aria2`。这些临时文件不能当成下载完成的文件；详见 `DATA_DOWNLOADS_zh.md`。

## 从真实数据建立试验集

```bash
python scripts/build_pilot.py --limit 1000
cat reports/pilot_summary.json

# 独立导入完整箭头数据；不尝试推断基元步骤分组。
mechbridge ingest mech-csv data/raw/mech_uspto/mech-USPTO-31k.csv \
  data/processed/mech_uspto_all.jsonl \
  --reaction-column updated_reaction --arrows-column mechanistic_label

# 导入 RGD1；按事件 ID 连接原子映射 CSV，先核对端点图一致性。
mechbridge ingest rgd1 data/raw/rgd1_zenodo/RGD1_allrxns.h5 \
  data/processed/rgd1_all.jsonl \
  --mapped-csv data/raw/rgd1_zenodo/RGD1CHNO_AMsmiles.csv

mechbridge pair data/processed/flower.jsonl data/processed/rgd1.jsonl \
  data/processed/pair_candidates.jsonl
```

试验集按源文件顺序选取，用来检验工程流程，**不是随机抽样，也不是泛化评测集**。精确重合可能为零；这不等于桥接研究不可行。没有重合的事件仍可通过补算电子结构和人工标注建立配对。

## 运行三维计算

```bash
# 实际的量子化学管线测试：构造的 NH3 反转，HF/STO-3G。
python scripts/qc_smoke.py

# 从下载的 Transition1x 数据导出一个真实反应事件。
mechbridge export-xyz data/processed/transition1x.jsonl \
  transition1x:C2H2N2O/rxn2091 work/rxn2091

# 明确声明电荷和电子态，再计算这个 TS 的电子结构描述。
mechbridge ibo work/rxn2091/ts.xyz work/rxn2091/ibo \
  --charge 0 --multiplicity 1 --method wb97x --basis '6-31g(d)'

# 高成本步骤：在所指定的 PES 上重新搜索并检查。
mechbridge verify work/rxn2091/ts.xyz work/rxn2091/verification \
  --dimer --charge 0 --multiplicity 1 --method wb97x --basis '6-31g(d)'

# 准备 10 条小体系计算任务；只生成任务，不批量启动昂贵计算。
python scripts/prepare_jobs.py data/processed/transition1x.jsonl work/jobs \
  --limit 10 --max-atoms 12 --charge 0 --multiplicity 1
```

`verify` 的内置双侧过程是**负模态位移 + 极小化，不是严格 IRC**。它会报告驻点力、投影频率、双侧极小化状态和电子能垒；`expected_endpoint_match` 和 `distinct_basins` 保持待核实。距离矩阵不能区分镜像等情形，不足以判定机理。没有找到预期 TS 不会写成“不可能反应”。

要接 MLIP，实现 `module:function` 工厂，接收 `charge`、`multiplicity` 及模型配置，返回新的 ASE Calculator：

```bash
mechbridge verify seed.xyz work/mlip --charge 0 --multiplicity 1 \
  --calculator my_backend:create_calculator --calculator-config model.json --dimer
```

工厂必须核对元素、电荷、自旋和模型适用域；本项目不随意中和体系。当前未下载、安装或运行 AIMNet2-rxn、FlowER 模型权重；4060 上也不要假定能直接运行 FlowER 官方完整模型。

## 数据与物理解释边界

- mech-USPTO 是模板生成的原始标签。`updated_reaction` 含上游补入试剂；项目只是保留上游记录，没有再次运行 MechFinder。原始行和标签来源均保留。
- 箭头序列不等于基元事件分组。形如 `301.1` 的源特定小数位点保留为未解析引用，不截断成原子 301，不擅自生成氢原子。
- RGD1 的 `RG/PG` 在作者解析说明中为未优化几何；所列 R/P 能量可能对应单独优化的碎片之和，不能把能量和该坐标冒充同一次单点计算。
- RGD1 映射 CSV 的图与 HDF5 的端点图先核对；映射编号与三维原子顺序仍单独待验证。
- Transition1x 的原始中间构型包含 NEB 优化采样；原始反应物、TS、产物三个结构也不能直接叫一条 IRC。
- 普通 MLIP 输出不包含轨道。IBO 分析需要波函数；轨道连续性存在表示和置换不确定性。当前没有硬规则箭头解码器。
- `compute_ibo` 默认各帧独立。只有输入真实有序路径时才加 `--ordered-path`；R/TS/P 三个离散构型不应假称细采样路径。低重叠会标记，不强行声称唯一轨道跟踪。
- 源标签、重新计算、模型输出分开保存；当前气相闭壳层代码不覆盖自由基、多重态交叉、光化学、金属催化和溶液标准态。
- 源方法名称一致不代表不同程序、网格、设置完全一致；PySCF 补算不能自动升级为源 ORCA/Gaussian 计算的完全复现。
- 当前 split 防止同端点对和同 sequence 泄漏，但尚未实现所有 scaffold、反应家族及跨源相同网络的完整隔离。

## 下一步真正需要补的研究工作

1. 从物理池筛选一批反应，补同层级鞍点优化、严格 IRC 与端点身份确认。
2. 沿已验证的路径补电子结构，人工审核首批箭头与轨道变化的对应；保留等价解释。
3. 定义源/目标电子位点的完整序列模型、三维种子生成器和反向解释模型，再用少量高可信配对数据训练。
4. 将计算成功、意外有效事件、未解决和超适用域分别作为反馈；比较相同计算预算下的有效事件数。
5. 再加入溶剂、标准态、构象与多分子微观动力学。

相关来源与许可见 `docs/DATASETS.md`；实际运行结果见 `reports/`。本项目不会把“框架跑通”写成“有机反应机理预测已经完成”。
