# 核心数据下载（本机版）

AIMNetCentral 广元素评测已下载 `aimnet2`、`aimnet2-2025`、`aimnet2-nse` 和当前注册表 `aimnet2-rxn` 的 member0，每个约 8.8–9.0 MB，保存于 `models/aimnetcentral/`。来源 URL、字节数与 SHA256 见 [model_receipt.json](reports/aimnetcentral_v8/model_receipt.json)。这些权重不计入下方数据集文件清单。

双分子验证新增 **Coley 组 [3+2] 环加成数据**：[作者仓库](https://github.com/coleygroup/dipolar_cycloaddition_dataset)，固定提交 `1608a64499779bf8f89a886308d97c6a830d084c`。已下载数据表（2,185,526 字节）和参考结构档案（94,023,449 字节），校验见 [下载回执](data/raw/coley_dipolar/receipt.json)。复现命令：`python scripts/data/fetch_coley_dipolar.py --reference-profiles`。3217、3216 的参考 TS 仅用于单独诊断，反应物起始主实验未读取参考 TS/产物几何。这些文件不计入下方原有五文件清单。

网络探索新增资源：SynEPD v0.4.1 符号数据和 AIMNet2-rxn 预训练权重已另存于 `data/raw/synepd/`、`models/aimnet2-rxn/`。固定提交与 SHA256 见 [探索资源清单](reports/exploration_resource_receipt.json)，复现下载用 `python scripts/data/fetch_exploration_resources.py`。这些新增文件不计入下方原有五文件清单。

**本机实际完成状态见 [local_download_inventory.json](reports/local_download_inventory.json)。只有 `all_verified: true` 才表示五个文件全部完成；旧的 `reports/data_inventory.json` 提供预期大小和哈希，其历史 `status` 不代表本机状态。**

| 数据              | 作者公开下载地址                                                                      |    原始字节数 | 项目内保存位置                                 |
| ----------------- | ------------------------------------------------------------------------------------- | ------------: | ---------------------------------------------- |
| mech-USPTO-31K v2 | [CSV](https://ndownloader.figshare.com/files/44708185)                                 |    32,348,446 | `data/raw/mech_uspto/mech-USPTO-31k.csv`     |
| FlowER v3         | [data.zip](https://ndownloader.figshare.com/files/55904909)                            |   238,206,912 | `data/raw/flower/data.zip`                   |
| RGD1 作者镜像     | [映射 CSV](https://zenodo.org/api/records/7860446/files/RGD1CHNO_AMsmiles.csv/content) |    60,724,906 | `data/raw/rgd1_zenodo/RGD1CHNO_AMsmiles.csv` |
| RGD1 作者镜像     | [HDF5](https://zenodo.org/api/records/7860446/files/RGD1_allrxns.h5/content)           | 1,340,003,744 | `data/raw/rgd1_zenodo/RGD1_allrxns.h5`       |
| Transition1x v4   | [HDF5](https://ndownloader.figshare.com/files/36035789)                                | 6,623,929,000 | `data/raw/transition1x/Transition1x.h5`      |

合计 **8,295,213,008 字节（8.30 GB，7.73 GiB）**。直接获取完整 HDF5，无需此前的 gzip 分卷或合并清单。FlowER 数据以 ZIP 保存，原项目支持直接读取其中的文本，不必展开全部约 5.97 GB 内容。

在项目根目录用 PowerShell 执行：

```powershell
python scripts/data/download_core.py --connections 4 --workers 3
```

下载器使用项目内 `tools/aria2/aria2c.exe`，来自 [aria2 官方 1.37.0 Windows 64 位发布包](https://github.com/aria2/aria2/releases/tag/release-1.37.0)，许可位于 `tools/aria2/COPYING`。其他平台安装 aria2 后通过 `--aria2 /path/to/aria2c` 指定可执行文件。

每个文件只有通过 **大小、MD5、SHA256** 三项核对后才去掉 `.part` 后缀，生成同目录的 `.receipt.json`。中断后运行同一条命令继续；保留 `.part` 和 `.part.aria2`，并行下载中的文件逻辑大小不能用来判断完成百分比。已有下载进程时不要重复启动，也不要用原单连接脚本续传 aria2 的分块文件。

查看状态和每个文件的速度：

```powershell
Get-Content reports/local_download_inventory.json -Encoding UTF8
Select-String -Path reports/download_Transition1x.h5.log -Pattern 'DL:' | Select-Object -Last 1
Select-String -Path reports/download_RGD1_allrxns.h5.log -Pattern 'DL:' | Select-Object -Last 1
```

下载结束后可独立重验；此命令会重新读取全部 8.30 GB 文件：

```powershell
python scripts/data/download_core.py --verify-only
```

`.part` 文件、下载日志以及旧报告都不能证明下载成功。下载器使用 IPv4，低速连接由 aria2 重试；Figshare 的短时签名失效时重新请求公开入口。只要文件仍有写入进展便继续尝试，连续五次失败且没有写入进展则记录失败。后台任务 PID 记录在 `reports/download_process.json`。

数据版本与许可见 [DATASETS.md](docs/DATASETS.md)。代码评估、文献及模型选择见 [FEASIBILITY_AND_MODELS_zh.md](docs/FEASIBILITY_AND_MODELS_zh.md)。本次核心数据任务不包含模型权重。
