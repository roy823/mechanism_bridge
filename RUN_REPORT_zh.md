# 本次运行报告

> 本文件是随项目附带的历史运行报告。2026-09-29 本机重新下载的实际状态见 `reports/local_download_inventory.json`；本机代码检查见 `reports/local_code_audit.json`，分析见 `docs/FEASIBILITY_AND_MODELS_zh.md`。以下量化计算与数据获取结果不代表本机已经重新执行。

日期：2026-09-29（Asia/Shanghai）

## 已完成

- 4 个核心数据集的 5 个原始文件，共 8,295,213,008 字节，大小和 MD5 校验通过，SHA256 已记录。
- mech-USPTO-31K：读取到 31,364 行、63 个机理类别。
- RGD1 HDF5：读取到 176,992 个记录组；镜像 HDF5 与 Figshare v6 的 MD5 相同。
- 四库分别导入 1,000 条真实记录，共 4,000 条。mech-USPTO 试验集箭头解析失败数：0。
- FlowER 试验集的组成、电荷及完整原子映射守恒检查均为 1,000/1,000。
- RGD1 CSV 与 HDF5 允许方向反转后的严格端点图核对：{'False': 443, 'True': 557}。不一致记录保留隔离字段，不强行赋予对应。
- 13 项自动测试通过，包括刚体模态投影、Hessian 正负号、反向数据配对、拆分泄漏以及电子能垒禁止作为自由能输入。
- ASE Dimer 的 HF/STO-3G NH3 反转数值测试：优化收敛=True；检查状态=index1_with_two_relaxed_minima。
- 构造 NH3 测试：TS 负频数=1；TS 最负频率约 -1081.3 cm⁻¹；双侧均收敛至无显著虚频结构。该构造例只验证程序流程。
- 对下载的真实事件 `transition1x:C2H2N2O/rxn2091` 的 TS，执行 PySCF ωB97X/6-31G(d) SCF 与 IBO/IAO 分析：SCF 收敛=True，轨道局域化收敛=True。计算得到的是电子结构特征，尚无人工箭头标注。
- 流式能量/力加载器通过真实 Transition1x 首帧检查。
- 准备 10 个小体系计算任务，未自动批量提交昂贵计算。

## 尚未完成

- 当前试验集跨库严格端点匹配候选：0。这仅是前 1,000 条的检索结果，不能外推为整个数据库没有重合。
- 经完整连接性和箭头证据验证的配对事件：0。
- 未执行所有有机 TS 的重新优化、Hessian 和严格双侧 IRC；当前双侧后端为负模态位移+BFGS。
- 未实现端点三维原子顺序与符号映射的自动认证、等价箭头集合、训练后的反向解码器或多分子溶液动力学。
- 未训练生成模型，未运行 FlowER/AIMNet2-rxn 权重。数据来源不代表模型已安装。
- PMechDB/PMechRP 未下载：官方入口需要填写个人信息及接受 CC-BY-NC-ND 条款。

## 可复现证据

原始文件清单：`data_inventory.json`；处理结果：`../data/processed/`；数值测试日志：`tests.log`；量化计算：`qc_smoke/`、`dimer_smoke/`；真实 TS 的电子结构：`real_event/`；任务列表：`prepared_jobs/jobs.jsonl`。

所有配对候选仍为 `verified_pair=false`。不应把本报告描述成已经完成机理真值数据库或训练成功的预测模型。

## 完整数据交付校验

RGD1 三维库以 gzip 无损压缩提供；Transition1x 以 15 个 gzip 分卷提供。已实际运行合并脚本，恢复得到 6,623,929,000 字节，整文件 SHA256 与原始下载文件完全一致。详细结果见 `reports/restore_verification.json`。完整文件下载见单独的 `DATA_DOWNLOADS_zh.md`。
