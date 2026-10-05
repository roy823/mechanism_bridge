# JCTC 工作进展（分支 `jctc-progress`，2026-10-05）

本目录汇总 `claude/jctc-code-fixes` 上的代码修改与主战役结果，供合作者同步。
详细数字与作业记录见 `results_log.md`；文中出现的 `docs/jctc/results/...`、`$W/...` 是工作副本和集群上的路径，本分支只附了主要图（`figures/`）。

## 1. 在做什么

从一个反应物出发，在机器学习势 AIMNet2-rxn 上自动构建反应网络。循环如下：
1. 用化学规则提出可能的反应：已发表的箭头推电子模板（SynEPD）、14 条人工规则，以及共振式；
2. 构造初始几何，做 TS 搜索（Dimer → Sella P-RFO）；
3. 检查 Hessian；
4. 双向 IRC，确认两端都是极小；
5. 登记节点与边。

每个运行的预算是 16000 次 MLIP 计算。协议冻结为 FP-JCTC-1（`configs/FP-JCTC-1.json`，提交 9734566）。所有判据在运行前写进 `campaign_plan.md`，事后分析单独标注。

四种策略只在"构造初始几何"一步不同：
- geometry：随机扰动；
- center_random：只知道反应中心原子；
- bond_edits：用提案的净变键；
- arrows：净变键加箭头细节（孤对电子方向、推拉顺序、相遇取向）。

注意 bond_edits 与 arrows 用的是同一套箭头模板产生的提案。此外有两类对照：
- 给定反应中心（oracle）：衡量搜索流程本身的上限；
- pyGSM（单端 GSM）：同一 MLIP、同一预算的外部基线。

## 2. 主要结果

| 图 | 内容 | 结果 |
|---|---|---|
| 图 3（`figures/fig3_main.png`） | 30 个体系的探索效率 | 箭头模板提案让效率提高约 3 倍：bond_edits 3.08、arrows 2.55 个与根相连的物种对，geometry 0.88；但 arrows 低于 bond_edits（p = 0.013） |
| 图 5a/5c（`figures/fig5c_main.png`） | 箭头种子特征 | 历史种子好于新版（p = 0.004）；箭头反向或打乱后反而更好（Holm p = 0.039 / 0.023） |
| 图 5b（`figures/fig5b_pilot_selected.png`） | 同一净变键、不同箭头 | 酮 → 烯醇：已发表模板的写法有 47% 到达预测产物，bond_edits 为 7%；Coley 组没有信息（`fig5b_coley_overview.png`） |
| 图 4（`figures/fig4_final.png`） | 找回参考反应 | T1x 盲搜 arrows 13.4%、bond_edits 9.1%；**给定反应中心：本流程 73%，pyGSM 55%**（Coley 上 32% 对 10%）；符号库覆盖率 T1x 6%、Coley 86%、RGD1 3% |
| 图 6 | 随机 80 条边的 DFT 复核（ωB97X/6-31G(d)） | 23/80 确认（29%，下界：31 条卡在 DFT TS 优化 100 步上限）；确认的边势垒 MAE 约 4.8 kcal/mol |

## 3. 诊断：箭头种子为什么落后（事后分析，见 `results_log.md` §11）

- 代码里没有符号错误。箭头信息在几何上是对的：T1x 上箭头种子更接近参考 TS（10/13，p = 0.021）。
- 问题集中在"孤对电子方向"这一项约束。同一提案、同一随机数下：
  - 这一项生效时，到达预测产物的次数是 arrows 290 对 bond_edits 387（p = 1×10⁻⁷）；
  - 不生效时两者打平，箭头还略好。
  - 箭头反向或打乱之所以"更好"，是因为恰好关掉了这一项。
- 结论：箭头作为提案来源有用；FP-JCTC-1 里孤对电子方向项的实现方式有害。

## 4. 下一步

1. FP-JCTC-2：
   - 修正或去掉孤对电子方向项，事先登记后在留出数据上验证；
   - 同时加一个不含箭头的提案对照（同一引擎，用随机变键枚举）。
2. TransitionNet：
   - 用符号种子（箭头 + 变键）高效建网；
   - 再用一个新模型判断给定条件（温度、溶剂、pH）下哪些边可达（参考 ReactionAtlas，arXiv 2606.30778）；
   - 以 ReactionAtlas 公开数据作对标；
   - 目标体系是从小分子出发的糖酵解网络（含磷酸与溶剂）。
3. 潜空间动力学：把网络中的极小、TS 与 IRC 路径放进潜在表示空间，把反应理解为能量曲面上的轨迹（正在调研）。

## 5. 文件

- `results_log.md`：全部结果与诊断（§1–§12）；
- `campaign_plan.md`：战役设计与事先写定的判据；
- `methods_draft_FP-JCTC-1.md`：方法初稿（英文正文）；
- `pre_freeze_tests.md`：冻结前测试记录；
- `figures/`：主要图。
