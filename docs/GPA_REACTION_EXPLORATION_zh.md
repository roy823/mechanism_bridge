# 符号与几何引导的 GPA 式反应网络探索

> 本文记录 v1 协议。当前代码已进入位移幅度统一的 v2，见 [新验证方案](VALIDATION_V2_zh.md)。严格复现 v1 请使用 `reports/repository_snapshots/pre_validation_v2.zip` 中的源码；不要将当前代码的输出与 v1 当成同一协议。

## 当前实现

这是一套已经可以运行的分子反应探索器，沿用 GPA 的组织方式：从已知极小值生成种子，寻找一阶鞍点，沿负模向两侧下降，登记实际极小值和连接，再从根节点连通分量里的新节点继续探索。没有移植晶体空隙、周期性离子占据或固定宿主的几何假设。

搜索输入没有参考 TS、产物三维构型或能垒。符号库的箭头可回放得到一个假设产品图，用于提议，但端点不被强制匹配该图。MLIP 找到的意外有效连接保留；由 A 提出的种子如果实际连接 B/C，只登记 B/C，不能伪造 A/B 边。

## 直接复用的数据

- RGD1：已下载的 226 条小体系候选池，只导出反应物坐标、元素及明确的闭壳层假设。原池曾按源 R/P 图一致性等条件筛选，因此本次样本有既有筛选偏差，不代表所有 RGD1。
- SynEPD：作者 GitHub 固定提交 `4fefc016fd4d4e4305dec92586e593bafe3b3e5b` 的 `polar.json`、release manifest 和软件许可。manifest 给出 v0.4.1、数据 CC BY 4.0；1,926 个反应、8,129 根箭头。与先前 v0.4.0 网页 8,123 根的差异属于版本差异。
- 采用 CHNO、净中性、反应物单片段、两电子箭头范围；箭头在反应物上回放并核对公开产品图，通过后建立正反两个方向的提议库。
- 1,926 条源记录中，1,850 条在上述起步范围外，67 条通过回放，8 条解析/回放不支持，1 条没有不同的化学图。被本实现拒绝不代表源数据错误。
- 在既有 226 条 RGD1 小体系记录中匹配到三个起点：硝基甲烷 `MR_460039_0`、乙烯醇 `MR_4792_0`、乙醛 `MR_8342_1`。后两者属于同一化学家族。只匹配反应物，无需两个库的反应产物相同。

记录见 `reports/network_data_audit.json`；下载版本与 SHA256 见 `reports/exploration_resource_receipt.json`。

## 使用的模型与计算

- 预训练模型：AIMNet2-rxn，HF revision `13e9de20a1e878a18a371e94dfb2176629281cb6`，`ensemble_0.safetensors`。没有训练新模型，没有使用随机森林。
- 推理：官方 `aimnet==0.2.0`，PyTorch 2.8 CPU，单个成员，2 个 PyTorch 线程。模型前向和依赖中的邻居表辅助函数明确采用 eager 模式。
- 使用固定权重配置声明的长程库仑项；按模型卡 `needs_dispersion=false` 显式关闭额外 D3。0.2.0 的 family defaults 会覆盖该项，故通过官方参数显式指定，避免重复色散。
- ASE Dimer：种子位置与方向均来自起点和提议；Dimer separation 0.005 Å、最大平移 0.1 Å、最多 100 步。
- 验证：最大原子力 0.03 eV/Å；有限差分 Hessian 步长 0.005 Å，投影整体平移和转动；虚频计数阈值 -30 cm⁻¹，要求恰好一个显著虚频。
- 连接：负模正反位移 0.15 Å（全坐标向量范数），无约束 BFGS 最多 150 步，极小值同样检查力和 Hessian，要求没有显著虚频，能量不高于鞍点。
- 上述证据标记为 `MLIP_two_sided_mode_displacement_descent`，`is_IRC=false`、`DFT_verified=false`。普通下降不能冒充严格 IRC；能垒是这一 MLIP 势面上的电子能量差。
- 只检查有限差分和旋转下的数值一致性，不因此声称模型与 DFT 一致。检查记录为 `reports/network_potential_check.json`。

## 三组比较

| 策略 | 可用信息 | 种子生成 |
|---|---|---|
| geometry | 当前构型 | 随机内部原子对距离扰动，去掉刚体平移和转动 |
| bond_edits | 当前构型 + 公开符号假设的净变键 | 基于共价半径的目标键长、未变键几何保持、排斥近距离碰撞 |
| arrows | 与上一组相同 + 电子源—汇配对 | 在同一个几何拟合中增加源键与汇键的进度同步残差 |

所有偏置只用于构造种子，最终 Dimer、Hessian 与两侧下降使用同一个未加偏置的 MLIP 势面。源—汇同步是可检验的启发式；异步协同反应可能不受益，不能当作电子转移的普适定律。当前箭头策略尚未利用孤对的真实空间方向或量子电子响应。

每个起点/策略使用同一个协议：最多 12 次搜索，每次最多 700 次势评估，总计最多 6,000 次（含初态优化及所有有限差分验证）。每批每节点 3 个种子；优先扩展未访问的新节点，随后在已覆盖节点使用新随机扰动继续搜索。所有失败计算也占预算。实际消耗另行报告，不把相同上限说成相同实际消耗。

节点按化学图及允许对称原子置换的 proper-rotation RMSD、能量差匹配，容差 0.15 Å 与 0.03 eV。TS 采用端点对、能量和固定原子身份下 RMSD 的保守去重；对称等价 TS 可能仍重复，化学图对数量另报。网络保留不同 TS 的平行边。

当前符号库仅支持完整体系精确图匹配。没有匹配箭头的新节点明确记为 unsupported；没有暗中调用参考产物补全规则。这是探索覆盖限制，不是反应不可行结论。

## 复现

Windows Python 3.11，推荐全新虚拟环境：

```powershell
python -m venv .venv-explore
.venv-explore\Scripts\python.exe -m pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu
.venv-explore\Scripts\python.exe -m pip install -r requirements-exploration-lock.txt
.venv-explore\Scripts\python.exe scripts/data/fetch_exploration_resources.py
.venv-explore\Scripts\python.exe scripts/data/prepare_network_pilot.py
.venv-explore\Scripts\python.exe scripts/exploration/check_reaction_potential.py
.venv-explore\Scripts\python.exe scripts/exploration/run_network_exploration.py --outdir reports/network_exploration_v1 --attempts 12 --seeds-per-node 3 --evaluations 6000
.venv-explore\Scripts\python.exe scripts/exploration/summarize_network_exploration.py reports/network_exploration_v1
```

需要既有 `data/processed/feasibility_pool.jsonl`，其构建脚本是 `scripts/data/prepare_feasibility.py`。运行器拒绝覆盖结果目录，复现实验请选新名字。源代码哈希、包版本、模型哈希、协议和随机种子随每次运行写入 manifest。

数据文件 `network_starts.jsonl` 只含起点，搜索器无权通过参数获得参考 TS/P。每个尝试保存 seed、初始方向、Dimer 轨迹、TS、下降轨迹及验证 JSON；网络可视化和表格从这些实际输出生成。

## 结果解释与后续

分子可视化使用 `scripts/exploration/render_molecular_exploration.py reports/network_exploration_v1` 生成，入口为 `reports/network_exploration_v1/molecules/index.html`。包含实际端点结构式、分子式、TS/端点三维球棍模型、逐帧下降序列、能量曲线、关键距离、结构式网络，以及 XYZ/MOL 和 PNG/PDF 导出。只使用已保存的构型和能量，未重新计算或插值；显示时整条序列采用同一个刚体变换。三维库在本地保存，页面无需联网加载资源。

看 `reports/network_exploration_v1/RESULTS_zh.md`、`summary.json` 和 `index.html`。重点分别报告根节点连通边、化学图对、新物种与从新节点继续探索，而非把所有构象边当作新化学反应。

本轮是三个起点、两个化学家族的工程与方法可行性试验。AIMNet2-rxn 训练包含 RGD1，不能当作未见化学泛化测试。反向符号假设来自公开箭头回放，没有从新波函数反推箭头，也没有证明新的物理标签改善学习模型。

后续应扩大反应家族与随机重复，补全可迁移的局部符号提议，抽样做 DFT/IRC，再比较同等预算下净变键和完整箭头的额外价值。新增计算证据可以回写事件记录，当前不自动改写外部符号数据库。
