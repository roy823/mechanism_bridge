# 核心数据下载（本机版）

原文的 `sandbox:/workspace/...` 是其他环境的临时附件地址，不能用于本机下载。本页改用作者公开的固定版本文件。

**本机实际完成状态见 [local_download_inventory.json](reports/local_download_inventory.json)。只有 `all_verified: true` 才表示五个文件全部完成；旧的 `reports/data_inventory.json` 提供预期大小和哈希，其历史 `status` 不代表本机状态。**

| 数据 | 作者公开下载地址 | 原始字节数 | 项目内保存位置 |
|---|---|---:|---|
| mech-USPTO-31K v2 | [CSV](https://ndownloader.figshare.com/files/44708185) | 32,348,446 | `data/raw/mech_uspto/mech-USPTO-31k.csv` |
| FlowER v3 | [data.zip](https://ndownloader.figshare.com/files/55904909) | 238,206,912 | `data/raw/flower/data.zip` |
| RGD1 作者镜像 | [映射 CSV](https://zenodo.org/api/records/7860446/files/RGD1CHNO_AMsmiles.csv/content) | 60,724,906 | `data/raw/rgd1_zenodo/RGD1CHNO_AMsmiles.csv` |
| RGD1 作者镜像 | [HDF5](https://zenodo.org/api/records/7860446/files/RGD1_allrxns.h5/content) | 1,340,003,744 | `data/raw/rgd1_zenodo/RGD1_allrxns.h5` |
| Transition1x v4 | [HDF5](https://ndownloader.figshare.com/files/36035789) | 6,623,929,000 | `data/raw/transition1x/Transition1x.h5` |

合计 **8,295,213,008 字节（8.30 GB，7.73 GiB）**。直接获取完整 HDF5，无需此前的 gzip 分卷或合并清单。FlowER 数据以 ZIP 保存，原项目支持直接读取其中的文本，不必展开全部约 5.97 GB 内容。

在项目根目录用 PowerShell 执行：

```powershell
python scripts/download_core.py --connections 16 --workers 3
```

下载器使用项目内 `tools/aria2/aria2c.exe`，来自 [aria2 官方 1.37.0 Windows 64 位发布包](https://github.com/aria2/aria2/releases/tag/release-1.37.0)，许可位于 `tools/aria2/COPYING`。其他平台安装 aria2 后通过 `--aria2 /path/to/aria2c` 指定可执行文件。

每个文件只有通过 **大小、MD5、SHA256** 三项核对后才去掉 `.part` 后缀，生成同目录的 `.receipt.json`。中断后运行同一条命令继续；保留 `.part` 和 `.part.aria2`，并行下载中的文件逻辑大小不能用来判断完成百分比。已有下载进程时不要重复启动，也不要用原单连接脚本续传 aria2 的分块文件。

查看状态和每个文件的速度：

```powershell
Get-Content reports/local_download_inventory.json -Encoding UTF8
Get-Content reports/download_Transition1x.h5.log -Tail 10 -Encoding UTF8
Get-Content reports/download_RGD1_allrxns.h5.log -Tail 10 -Encoding UTF8
```

下载结束后可独立重验；此命令会重新读取全部 8.30 GB 文件：

```powershell
python scripts/download_core.py --verify-only
```

`.part` 文件、下载日志以及旧报告都不能证明下载成功。Figshare 单文件接口会发放短时签名链接，下载器在失败后重新请求公开入口并利用 aria2 续传。

数据版本与许可见 [DATASETS.md](docs/DATASETS.md)。代码评估、文献及模型选择见 [FEASIBILITY_AND_MODELS_zh.md](docs/FEASIBILITY_AND_MODELS_zh.md)。本次核心数据任务不包含模型权重。
