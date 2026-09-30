# 项目可行性、文献与模型选择

## 2026-10-01：AIMNetCentral 广元素实测选择

当前官方维护入口为 [AIMNetCentral](https://github.com/isayevlab/aimnetcentral)。本机已经固定并实测四个 member0：

| 用途 | 选择 | 元素范围 | 本项目结论 |
|---|---|---|---|
| CHNO 历史复现 | 原固定 `aimnet2-rxn` | H/C/N/O | 保持旧结果可比，不追溯替换推理默认值 |
| 广元素闭壳层探索 | **`aimnet2-2025`** | H/B/C/N/O/F/Si/P/S/Cl/As/Se/Br/I | 当前首选；分子间作用改进，速度与其他成员相当 |
| 较高参考层级的广元素计算 | `aimnet2` | 同上 | wB97M-D3；三个 DFT-seeded CHNO 事件均严格保留 |
| 自由基/开壳层 | `aimnet2-nse` | 同上 | 支持 multiplicity；只在需要时使用，当前检查成本较高 |
| Pd 催化 | `aimnet2-pd` | 14 元素中的 As 换成 Pd | 独立分支，尚未下载或评测 |

本机 11 原子体系上，四个模型 CPU 单构型约 8–9 ms；RTX 4060 单构型约 11–14 ms，未体现优势。batch32 时 CPU 约 772–781 构型/s，GPU 约 2609–2782 构型/s。广元素模型显存峰值约 106–128 MiB。搜索中的逐步 Dimer 保持 CPU；批量 Hessian、NEB 图像与较大体系优先 GPU。

在丙酮互变、甲醛水合、甲醛二聚三个 DFT TS/负模辅助检查中，`aimnet2` 与 `aimnet2-2025` 均严格保留 3/3 端点对；NSE 保留 3/3 图对，但甲醛二聚有一个端点没有通过极小值验收。官方文档的 Cl⁻ + CH₃Cl 对称 SN2 示例在三个广元素模型上均通过本项目的一阶鞍点和双侧下降检查。这些都是参考辅助检查，不是自主发现或全面精度 benchmark。

当前 AIMNetCentral 注册表会给 `aimnet2-rxn` 加外部 D3，而本项目历史固定 HF artifact 按当时模型卡显式关闭 D3。在测试几何上，两种默认设置产生 0.2146 eV 能量差和 0.01836 eV/Å 最大力差；关闭注册表 D3 后两者完全一致。因此旧网络不静默升级。

完整数字和权重哈希见 [AIMNetCentral v8 报告](../reports/aimnetcentral_v8/RESULTS_zh.md)。

> 历史评估记录，保留当时的数据审计与判断。当前实现、模型安装状态和入口以 [README](../README_zh.md) 为准；最新试验见 [验证 v2](VALIDATION_V2_zh.md)。

评估日期：2026-09-29。依据为本机全部 `src/mechbridge` 模块、下载与试验脚本、处理样本和公开一手文献。本机下载状态以 `../reports/local_download_inventory.json` 为准。

用户补充 GPA 项目后，研究主线进一步明确为“几何/势能发现事件—电子结构—符号解释”。已有能力、迁移边界及两个模型的分工见 [GPA_BRIDGE_DESIGN_zh.md](GPA_BRIDGE_DESIGN_zh.md)。下面的跨库零交集仅是一种严格记录检索结果，不是项目可行性的前提，也不表示缺少共同机理类型。

**判断：适合继续做“小体系、明确适用域内的机理—三维路径证据配对”研究原型；目前尚不具备端到端机理预测、自动弯箭头真值构建或实验产率预测能力。** 首要工作是获得可审核的真实配对数据。单独拼接现成模型无法解决配对和验证缺口。

**代码已经完成什么**

| 模块 | 本次读代码确认的能力 | 后续用途 |
|---|---|---|
| `adapters.py`、`io.py` | 四库导入、来源与原始记录保留、坐标和单位转换 | 数据入口与溯源 |
| `chemistry.py` | 图标准化、组成/电荷/映射身份核对、净键变化 | 数据筛查；净键变化不唯一确定箭头 |
| `pairing.py` | 完整端点图检索、反向匹配；按端点和 sequence 联合划分 | 候选检索与部分泄漏控制 |
| `physics.py` | Dimer、差分 Hessian、刚体模态投影、两侧 BFGS | 局部驻点检查；连接性仍待补齐 |
| `electronic.py` | SCF、IAO/IBO、相邻帧轨道重叠匹配 | 电子结构特征；输出 `arrows=None` |
| `backends.py` | 气相闭壳层 PySCF；外部 ASE 工厂 | MLIP 接口存在，具体 AIMNet2 后端尚未实现 |
| `training_data.py` | 按源反应分组流式加载能量/力 | 数据加载器；未定义网络、损失和训练循环 |
| `kinetics.py` | 活化 Gibbs 自由能的单分子 TST 与速率矩阵 | 未覆盖溶液多分子网络、浓度与产率 |

本机执行 `python -m pytest -q`：**13 passed in 1.25s**。这仅覆盖现有单元测试。本次未重新运行 PySCF 量化计算，随包提供的历史量化日志不是本机新结果。

**影响可行性的具体缺口**

1. **符号和三维数据尚未对齐。** 本次重新读取现有 4,000 条示例：FlowER 的 1,000 条有端点但无显式 `arrow_pairs`；mech-USPTO 的 1,000 条有箭头，但都是 `symbolic_overall_reaction`，基元分组未知，会被配对函数排除。RGD1 有端点；Transition1x 的 1,000 条均无端点 SMILES，会被配对函数跳过。当前不能直接完成四库互配。
2. **Transition1x 未进入试验配对与自建划分。** `build_pilot.py` 将其写盘后没有加入 `physical`；实际配对和 `pilot_split.jsonl` 只覆盖 FlowER/RGD1。能量/力训练应保留官方 split，并按反应而非帧划分。
3. **原子对应尚未认证。** 已有 RGD1 前 1,000 条中，CSV/HDF5 端点图 557 条一致、443 条不一致；即便一致，也未证明 SMILES 映射号与 XYZ 原子顺序对应。不可把行号或原子索引当作映射号。对 443 个冲突的诊断发现：去掉立体/同位素信息后 230 条图一致，另 213 条仍不一致；这仅用于查错，不能据此放宽正式配对条件。新下载的完整映射 CSV 有 176,898 个唯一事件 ID，需与 HDF5 核对覆盖差异。详见 `../reports/local_data_diagnostics.json`。
4. **缺少严格 IRC 和端点认证。** `verify_saddle()` 明确返回 `is_IRC=False`，`expected_endpoint_match='not_checked'`，也不认证两侧为不同势阱。需要同一势能面上的 TS 优化、有效单虚频、双向 IRC、端点优化及包含立体信息的身份匹配。可采用已有 [PySisyphus IRC](https://pysisyphus.readthedocs.io/en/latest/irc.html) 或成熟量化程序的实现。
5. **物理条件缺失。** 四组现有样本的多重度均未知。端点图一致不能替代电荷、自旋、溶剂、质子转移辅助分子和构象的核对，补充假设必须有记录。
6. **轨道到箭头仍是研究任务。** 匈牙利匹配给出一种轨道对应；近简并、局域化多解和低重叠仍有歧义。应保留等价箭头解释，并人工审核，不能仅凭键级差或人口阈值宣布获得唯一机理。
7. **量化成本高。** 中心差分 Hessian 对 N 个原子需要 6N 次梯度调用；TS 与两侧极小值约为 18N 次，还不含优化和 IRC。20 原子体系三次 Hessian 约需 360 次梯度。MLIP 适合前期筛选，最终标注仍需独立量化证据。
8. **评估隔离不足。** 现有 split 尚未完成跨库结构家族隔离及预训练重叠排查。历史候选数 0 只是前 1,000 条结果；本次额外的完整检索结果如下。

**本次完整重合检索**

`scripts/diagnostics/audit_overlap.py` 扫描原版 `flower_dataset` 的 train/val/test 共 **1,622,935 行**，与 RGD1 映射 CSV 的 **176,898 条记录、139,693 个非恒等端点图对**比较。跳过 360,463 行相同 SMILES 端点，以及 1,259,214 行超出 RGD1 的 CHNO/10 重原子范围的数据，剩余 3,258 行进行图标准化；严格候选数为 **0**。结果见 `../reports/full_overlap_audit.json`，可用 `python scripts/diagnostics/audit_overlap.py` 复现。

这一结果只针对原版 FlowER 和该映射 CSV、当前 RDKit 标准化规则及完整体系的精确匹配；未扫描新数据变体，也不排除 HDF5 原始端点、局部反应中心或新计算体系的重合。部分 RGD1 源 SMILES 会触发 RDKit 立体方向冲突警告，日志予以保留，结果不是经过化学人工审核的全库定论。它足以说明：不能把“直接跨库取交集”当成获得大量正配对的默认方案，应优先从物理路径建立箭头标注。

**文献支持什么**

| 一手来源 | 对应环节 | 证据边界 |
|---|---|---|
| Chen et al., Scientific Data 2024，[mech-USPTO-31K](https://doi.org/10.1038/s41597-024-03709-y) | 专家模板扩展得到 31,364 条箭头标注 | 是模板生成的机理假设，不是逐条 TS/IRC 证明 |
| Joung et al., Nature 2025，[FlowER](https://doi.org/10.1038/s41586-025-09426-9) | 电子重分配与 flow matching 的符号反应生成 | 另需三维与动力学验证；生成流时间不是实际反应时间 |
| Zhao et al., Scientific Data 2023，[RGD1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10025260/) | CHNO 小体系反应/TS 的物理数据基础 | 源验证层级、构象和端点能量对应仍须保留并核查 |
| Schreiner et al., Scientific Data 2022，[Transition1x](https://doi.org/10.1038/s41597-022-01870-w) | 约 960 万个反应区域能量/力计算 | NEB 优化采样不等于有序 IRC；没有配套弯箭头真值 |
| Knizia & Klein, Angew. Chem. Int. Ed. 2015，[IBO 与电子流](https://doi.org/10.1002/anie.201410637) | 沿反应坐标的 IBO 变化解释电子对迁移 | 不等于已有跨反应族通用、唯一的箭头解码器 |
| Anstine et al., 2025，[AIMNet2-rxn 预印本](https://www.cambridge.org/engage/chemrxiv/article-details/685505c9c1cb1ecda0f701de) | 反应专用 MLIP 加速路径与能垒筛选 | 本次核到的是预印本版本；需要误差与适用域评估 |
| Tuo et al., Nature Communications 2026，[MolGEN](https://doi.org/10.1038/s41467-026-75654-w) | 产物/TS 的条件 flow matching 和网络生成 | 生成产物+TS 已有先例；需证明符号—轨道配对或反馈的额外收益 |

据此，建议把研究问题聚焦于可追溯的“符号步骤—已验证路径—电子结构变化”配对，以及它能否改善候选筛选效率和解释可信度。该定位是本次分析的建议，尚无实验结果证明性能提升，也不构成穷尽性的新颖性查新。

**需要下载哪些模型**

| 优先级 | 文件/工具 | 用途和作者入口 | 本机限制 |
|---|---|---|---|
| 最初无需权重 | PySCF + ASE/RDKit | 从现成 TS 开始做电子结构；[PySCF 安装](https://pyscf.org/user/install.html) | 本机缺 PySCF；官方不支持原生 Windows，使用 WSL2/Linux |
| 优先模型 | `isayevlab/aimnet2-rxn` 的 `config.json`、`ensemble_0.safetensors`，保留模型卡和 revision | 能量/力与反应路径预筛；[官方模型仓库](https://huggingface.co/isayevlab/aimnet2-rxn) | 先限 CHNO、净中性、闭壳层；8GB 小体系试运行尚未实测 |
| 不确定性实验再补 | `ensemble_1.safetensors` 至 `ensemble_3.safetensors` | 同一 AIMNet2-rxn 的四成员集成 | 可逐个加载；成员分歧不是严格误差界 |
| 符号生成阶段 | FlowER v3 `checkpoints.zip`，444,244,997 字节 | [官方权重](https://ndownloader.figshare.com/files/55904912)；[代码](https://github.com/FongMunHong/FlowER) | 官方至少 25GB 显存；本机 8GB 不满足默认复现配置，建议远程 GPU |
| 三维生成阶段可选 | MolGEN TS checkpoint；做产物生成时再补对应 checkpoint | [官方 Model weights 入口](https://github.com/tuoping/MolGEN#model-weights) | 尚无适配器；需核对预处理、原子顺序与显存；初期使用库中 TS 即可 |

本次数据下载只取四库五个原始文件，**没有下载或运行上述模型权重**。当前流程不需要通用大语言模型。PySCF、PySisyphus 是计算软件，不是神经网络 checkpoint。MLIP 的能量/力不提供 IBO 所需的占据轨道，电子结构阶段仍需量化波函数；当前也没有可直接下载并接入的“轨道→弯箭头”解码器，需要自行定义、标注和训练。

AIMNet2-rxn 的参考层级是 ωB97M-V/def2-TZVPP，能量采用按元素平移的标度，应比较同模型、相同组成下的能量差。不能把其绝对能量与 Transition1x/RGD1 直接混合训练。训练数据包含重新计算的 RGD1，因此 RGD1 不能未经排重就用作该模型的独立泛化测试。模型卡同时存在旧版限制及 `NEXT_VERSION` 描述，工厂必须自行检查元素、电荷、自旋与有限数值。[模型卡](https://huggingface.co/isayevlab/aimnet2-rxn)

FlowER v3 权重包 MD5 为 `2ab43d1956aafc12ead0b707b2b3c9db`，含不同模型/数据变体。复现 Nature 论文应使用官方指明的 `2.0.0` 标签并匹配具体 checkpoint；主分支已转为新数据工作流。不能将 `flower_dataset` 与 `flower_new_dataset` 合并视为独立样本。缩小 batch 以适配 8GB 需要实测，不能承诺成功。[官方版本说明](https://github.com/FongMunHong/FlowER)

本机：Python 3.11.9、RTX 4060 Laptop 8,188 MiB、PyTorch 2.5.1+cu121、RDKit 2023.9.5、NumPy 1.26.4、SciPy 1.15.3、h5py 3.16.0、ASE 3.25.0。RDKit 低于项目要求的 2024.3；测试通过不证明所有化学边界行为一致。当前 AIMNet 主分支要求 Python ≥3.11、PyTorch ≥2.8，建议在独立 WSL2/Linux 环境安装并锁定版本。[AIMNet 安装说明](https://github.com/isayevlab/aimnetcentral)

需要 FlowER 权重时，项目现有脚本支持：

```bash
python scripts/data/download_data.py --datasets flower --checkpoints --max-file-mb 500
```

在独立 AIMNet 环境中，官方接口可获取并加载默认成员。以下尚未在本机执行：

```python
from aimnet.calculators import AIMNet2Calculator
calc = AIMNet2Calculator("isayevlab/aimnet2-rxn", compile_model=False)
```

正式实验保存库版本、模型 revision、checkpoint SHA256 和计算参数。当前 `load_factory()` 不会代替后端完成这些核查。

**建议的实施顺序与验收条件**

1. 完成数据校验和格式读取，保留 FlowER 原始训练/验证/测试划分。
2. 从物理库挑约 20–50 个中性闭壳层 CHNO 小体系，完成同参考层级的 TS、频率、双侧 IRC、端点身份和原子对应。这是建议的工程起步规模，不是统计充分性保证。
3. 对真实有序路径做 IBO 跟踪，人工标注、交叉审核并保留等价解释，先获得几十条确证配对再扩大规模。
4. 接入 AIMNet2-rxn 预筛。区分预期反应成功、替代反应成功、未解决、域外；未收敛不能作为“反应不可能”的负标签。
5. 按需增加 FlowER 和一种三维种子模型。在同等 DFT 梯度预算下比较直接搜索、MLIP 预筛、生成种子和电子结构反馈，报告确证事件数、端点匹配率、能垒误差与人工箭头一致性。
6. 留出未用于预训练和标注设计的反应/骨架家族，再扩展溶剂、构象、自由能和多分子动力学。当前证据不支持全域机理或产率结论。
