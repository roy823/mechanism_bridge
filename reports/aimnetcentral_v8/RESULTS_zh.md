# AIMNetCentral 广元素模型评测 v8

结论：广元素闭壳层探索首选 `aimnet2-2025`；需要更高参考层级时用 `aimnet2`；自由基/开壳层才启用 `aimnet2-nse`；CHNO 历史结果继续固定原 `aimnet2-rxn`。

| 模型 | 元素数 | CPU 单点/ms | CPU batch32/s | GPU 单点/ms | GPU batch32/s | GPU峰值/MiB | CHNO严格通过 | SN2 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| AIMNet2 | 14 | 8.66 | 779 | 11.22 | 2782 | 106.5 | 3/3 | True |
| AIMNet2-2025 | 14 | 8.14 | 781 | 12.26 | 2653 | 118.3 | 3/3 | True |
| AIMNet2-NSE | 14 | 8.67 | 781 | 10.61 | 2609 | 127.8 | 2/3 | True |
| AIMNet2-rxn | 4 | 8.33 | 772 | 13.84 | 2693 | 76.1 | 3/3 | 不支持 Cl |

## 如何使用

- Dimer 的逐步单构型调用：对这些 11 原子体系 CPU 更快；继续用 CPU。
- 批量有限差分 Hessian、NEB 图像或更大分子：GPU 吞吐明显更高。
- 当前 11 原子体系上，GPU 批量有限差分 Hessian 比 GPU 原生二阶导数更快；CPU 原生 Hessian略快。当前工作流保留批量有限差分。
- 四成员 ensemble 约增加四倍计算量；搜索阶段使用 member0，只有筛选/不确定性阶段再计算 ensemble。

## 反应区检查

三个模型间共有的 DFT TS/负模辅助事件是丙酮互变、甲醛水合、甲醛二聚。`aimnet2` 与 `aimnet2-2025` 都严格保留 3/3；NSE 保留图对 3/3，但二聚端点未全部通过曲率检查。该检查使用了参考 TS 和负模，不是自主发现。

官方文档的 Cl⁻ + CH₃Cl 对称 SN2 例子在三个广元素模型上都通过严格检查；AIMNet2-2025 为 210 次评估、约 1.28 s CPU。

## 重要版本差异

当前 AIMNetCentral 注册表对 `aimnet2-rxn` 加入外部 D3；相对历史固定 HF 推理，测试几何的能量差为 0.2146 eV、最大力差为 0.0184 eV/Å。关闭 D3 后两者数值完全一致。旧结果不能静默换成注册表默认。

## 边界

- 14 元素范围：H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I；不含碱金属/碱土金属及除 Pd 外的过渡金属。
- 当前符号动作库仍以中性闭壳层 CHNO 为主。换势能扩大了物理后端，尚未自动扩大电子机理提议覆盖。
- GPU 结果使用现有 CUDA 环境做工程计时；正式部署应建立干净锁定环境。未评测 `torch.compile`。
- 不同模型参考理论层级不同，能量不能混在同一网络或直接拼接训练。

数据：[summary.json](summary.json)、[model_receipt.json](model_receipt.json)、[comparison.png](comparison.png)。
