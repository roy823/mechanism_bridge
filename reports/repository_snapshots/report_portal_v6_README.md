# 统一报告与 v6 证据包

文件：`report_portal_v6_evidence.zip`（73,206,127 字节）。

SHA256：`7ad2344861925a30299690f1461cbf2305d47b4d5a10b88ffc800e21ebd03e36`。

包含统一 HTML、实际分子坐标与下载文件、页面链接到的原始记录、完整 v6 实验、生成代码、测试与来源清单。历史全量原始轨迹仍以此前各轮证据包为准；本包不包含模型权重和大型原始数据库。

当前工作区直接打开 `reports/index.html`。在其他位置查看时，解压到一个新的空目录，打开其中的 `reports/index.html`；不需要运行服务器或安装模型。

例如，在仓库根目录执行：

```powershell
Expand-Archive -LiteralPath reports/repository_snapshots/report_portal_v6_evidence.zip -DestinationPath .cache/report_portal_v6
```

然后打开 `.cache/report_portal_v6/reports/index.html`。请使用新的目标目录，避免覆盖已有研究文件。

逐文件哈希见同名 `report_portal_v6_evidence.json`。旧 HTML 外观另保存在包内的 `pre_report_portal_html.zip`，原始数值与轨迹未因页面整合而改变。
