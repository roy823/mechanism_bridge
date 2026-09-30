# 数据来源、版本与许可

> 以下获取结果来自随包附带的历史清单。本机重新下载状态以 `../reports/local_download_inventory.json` 为准，下载入口见 `../DATA_DOWNLOADS_zh.md`。

检索与获取日期：2026-09-29（Asia/Shanghai）。所有最终文件均按公开元数据的大小和 MD5 校验，并计算了 SHA256。精确结果见 `reports/data_inventory.json`。

| 数据 | 固定来源 | 本次获取 | 许可与科学用途 |
|---|---|---|---|
| mech-USPTO-31K | [Figshare v2](https://doi.org/10.6084/m9.figshare.24797220.v2) | `mech-USPTO-31k.csv`，32,348,446 字节；31,364 行，63 个 mechanistic_class | CC BY 4.0；Shuan Chen、Ramil Babazade 等。模板生成的符号机理假设，不是逐步 DFT 真值 |
| FlowER | [Figshare v3](https://doi.org/10.6084/m9.figshare.28359407.v3) | `data.zip`，238,206,912 字节；包含 flower_dataset 和 flower_new_dataset | 条目标注 MIT；Joonyoung F. Joung、Mun Hong Fong 等。基元步骤的映射 SMILES/sequence ID；不保证显式唯一箭头 |
| RGD1 | [作者 Zenodo 条目](https://doi.org/10.5281/zenodo.7860446)；[Figshare v6](https://doi.org/10.6084/m9.figshare.21066901.v6) | `RGD1_allrxns.h5`，1,340,003,744 字节；176,992 个组。`RGD1CHNO_AMsmiles.csv`，60,724,906 字节 | CC BY 4.0；Qiyuan Zhao 等。HDF5 MD5 与 Figshare `RGD1_CHNO.h5` 完全一致；CSV 是 Zenodo 的独立文件，不能声称等于 v6 CSV |
| Transition1x | [Figshare v4](https://doi.org/10.6084/m9.figshare.19614657.v4) | `Transition1x.h5`，6,623,929,000 字节 | 条目标注 MIT；Mathias Schreiner 等。三维构型、能量、力及端点/TS；优化采样不等于 IRC |
| PMechDB/PMechRP | [官方下载页](https://deeprxn.ics.uci.edu/pmechdb/download) | **未下载**：需要填表并接受许可 | 下载页为 CC-BY-NC-ND；包括教科书路径资源。未提交个人信息或申请 |

核心原始文件合计 **8,295,213,008 字节**。FlowER ZIP 的旧、新变体不能直接拼接并当作互相独立的数据；本次试验使用 `data/flower_dataset/train.txt` 的前 1,000 行。压缩包元数据的总展开大小为 5,973,556,710 字节。

另为核查格式获取了 [MechFinder 官方仓库](https://github.com/snu-micc/MechFinder) 的原始 USPTO_33K CSV，固定提交 `a080f8fec82ecc08f6eac60a62d61034caabc034`。它不是 mech-USPTO-31K 标签库，未用于替代后者，亦未运行模板补全代码。

## 原始出处

- Chen et al., *A large-scale reaction dataset of mechanistic pathways of organic reactions*, Scientific Data (2024). [DOI](https://doi.org/10.1038/s41597-024-03709-y)
- Joung et al., *Electron flow matching for generative reaction mechanism prediction*, Nature (2025). [DOI](https://doi.org/10.1038/s41586-025-09426-9)；[官方代码](https://github.com/FongMunHong/FlowER)。仓库主分支已有后续更新，使用时需同时锁定代码提交和数据变体。
- Zhao et al., *Comprehensive exploration of graphically defined reaction spaces*, Scientific Data (2023). [DOI](https://doi.org/10.1038/s41597-023-02043-z)；[官方解析代码](https://github.com/zhaoqy1996/RGD1)
- Schreiner et al., *Transition1x—a dataset for building generalizable reactive machine learning potentials*, Scientific Data (2022). [DOI](https://doi.org/10.1038/s41597-022-01870-w)；[官方代码](https://gitlab.com/matschreiner/Transition1x)
- Knizia & Klein, *Electron Flow in Reaction Mechanisms—Revealed from First Principles*, Angew. Chem. Int. Ed. (2015). [DOI](https://doi.org/10.1002/anie.201410637)
- PySCF [IBO/IAO 官方实现](https://pyscf.org/_modules/pyscf/lo/ibo.html)

## 授权与再分发

本项目保留了原始数据条目的作者、DOI、许可声明、文件级校验值。原始数据仍受其各自许可约束，项目代码的打包不改变原始数据的权利归属。FlowER 上游整合来源应继续追溯，条目许可不能自动解决每个上游来源的所有使用条件。

## 验证等级

双分子扩展另使用 [Coley 组的 1,3-偶极环加成公开数据](https://github.com/coleygroup/dipolar_cycloaddition_dataset)。数据表和参考结构档案固定到提交 `1608a64499779bf8f89a886308d97c6a830d084c`，位于 `data/raw/coley_dipolar/`。主实验仅输入反应物；参考产物用于事后比较，参考 TS 的独立诊断不计入主搜索命中。完整说明见 [双分子验证](BIMOLECULAR_V5_zh.md)。

“下载成功”只表示字节校验成功。“导入成功”只表示该格式被解析。“图匹配”只是跨库候选。以上三者均不等于箭头—TS 配对已经得到量化计算验证。
