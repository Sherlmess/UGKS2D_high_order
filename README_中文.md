# 二维 S2 可复现底座：宏微观高阶重构＋二阶时间推进

发布版本：**`v1.0`**。数值底座身份：`S2D-O3-S2Steady-I3R1-20260923`。

[版本说明与后续计划](docs/RELEASE_v1.0_中文.md)

当前二维高阶专用仓库：**https://github.com/Sherlmess/UGKS2D_high_order**，SSH：`git@github.com:Sherlmess/UGKS2D_high_order.git`。原生隐式基础设施来自 [implicitUGKS](https://github.com/Sherlmess/implicitUGKS) 的 `v2.0`，提交 `753607a44572342afd094404ec545c005ab737a6`。本仓库独立维护整理后的二维 S2 代码；其来源与旧研究记录保留在 `provenance/`。

这是当前可复现的二维 **S2** 路线：宏观 W 与微观 h/b 均使用 WENO5-AO 重构，真实 UGKS 特征线面通量，低成本 S2 体碰撞，宏观显式更新＋微观自洽 CN 二阶时间推进；另含求同一离散稳态的 `S2Steady` 隐式加速。它不是旧的“宏观高阶、微观 Van Leer”路线，也不是 W+f＋2S4O 四阶时间算法。

## 从哪里开始

- [算法阅读与实现指南](docs/ALGORITHM_GUIDE_中文.md)：阅读顺序、公式和数据流。
- [配置、运行和输出](docs/RUN_GUIDE_中文.md)：7 个独立例子，不依赖原研究目录。
- [成果、限制及下一阶段](docs/RESULTS_AND_SCOPE_中文.md)：区分历史科学验证与本次整理检查。
- [冻结源码身份](provenance/SOURCE_IDENTITY.json)：32 个源码／构建文件逐文件哈希及来源。
- [此次验证结果](verification/RESULTS.json)：干净编译、控制流与保存态复现结果。

## 快速使用

需要 Python 3（脚本仅标准库）、GNU Fortran、Open MPI、make；macOS 还需可用的 Xcode Command Line Tools SDK。

```sh
git clone git@github.com:Sherlmess/UGKS2D_high_order.git
cd UGKS2D_high_order
git checkout v1.0
python3 -B scripts/verify.py
python3 -B scripts/build.py
python3 -B scripts/run_case.py bgk_saved_step --output runs/bgk_check_01
python3 -B scripts/run_case.py shakhov_saved_step --output runs/shakhov_check_01
python3 -B scripts/run_case.py steady_finish --output runs/steady_finish_01
```

每次选择新的输出目录，脚本拒绝覆盖旧运行。打包的源代码可从零编译；不携带依赖本机路径的可执行文件。`examples/` 的输入和预期二进制状态随包提供，旧研究目录路径只用作来源说明。

**当前实际限制**：单 MPI rank；均匀 `[0,2π]²` 双周期；N=8/16/32/64；GH16×16（256 个速度节点）；CK=1；源码固定 τ=0.1。Pr=1 是 BGK，已有 Pr=2/3 Shakhov 非定常证据。S2Steady 的完整收敛与成本资格限于 N8 BGK。不能直接把现有 S2 配置改成壁面、非均匀网格、多 rank 或变 τ 工况。

## 文件布局与维护

```text
src/             完整原生求解器源码：I3-R1 已验收版本，不需叠加历史补丁
Makefile         已验收构建文件，保持原字节
examples/        配置、网格、必要重启态和对应预期终态
scripts/         构建、隔离运行、只读校验
tests/           v2.0 原始跟踪测试＋本次终止逻辑单测
docs/            中文算法、运行、证据范围
provenance/      来源、历史结果摘要与数据
verification/    本次整理后的实际验证记录
SHA256SUMS       固定交付文件清单，排除本地产物 bin/obj/runs
```

旧 `tests/` 从不可变 v2.0 提取，并非父树未提交内容；本次没有重新执行完整 legacy 测试矩阵。原 Makefile 中这些目标保留兼容，但不能把其存在当成 S2 的额外科学验证。

此版本从独立完整快照建立，没有修改 `implicitUGKS_s2d1/src` 或旧研究证据。以后二维 S2 工作以本仓库为入口，不再从旧研究包逐层拼接补丁。`provenance/` 和 `verification/` 中的旧绝对路径、未提交字段是打包验证时的历史记录，不代表运行依赖或当前发布状态。
