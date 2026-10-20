# v1.1a 精简版：删除无调用代码，保留现有数值实现

来源为远端 `Sherlmess/UGKS2D_high_order` 的 `v1.1`，提交 `b8eb22b487fc6b4c42633f049dec089b4c7082c7`。本地另一个同名标签不是此次基线。此目录是可独立构建的精简版本；发布不改写 v1.1 历史，未提交超算作业。

## 修改与保留

源码由 22,940 行减少为 21,575 行，删除 1,365 行（5.95%）：删除 16 个没有仓内调用方的过程，以及 3 个失去消费者的模块，同步移除其构建依赖和导出声明。未新增数值公式、近似、参数或约束。

- 删除旧 `ConstantBGKModel2D` 及其点目标、单元平均目标例程；移除其专用的 `mod_CorrectedWF2DTypes`、`mod_WENO5AO2D`、`mod_WENO5AOScalar` 三个模块。
- 删除无调用的 `projector_direction_2d`、`stationary_e_2d`、`stationary_e_direction_2d` 以及随后失去消费者的 Hessian/方向导数辅助链。
- 删除 `target_at`、`s2d1_weno1d_companion`、`write_microVDF`、`bc_periodic` 和两个旧 `MPIUpdata_MacroVars_*dirc` 过程。

当前 S2 的生产 WENO5-AO 实现在 `mod_S2D1Reconstruction.f90` 中保留。UGKS 特征线面通量、S2 体求积/热流闭合、物理变 tau、PR(2)、微观前后单扫、质量比例校正、接受/拒绝门、MPI 接线、完整重启和仍有调用的诊断/历史入口全部保留。

36 个其余源码文件逐字节未变；7 个编辑文件中保留的 119 个过程体逐字节未变。修改仅涉及完整过程删除、旧类型删除、导出/导入及 Makefile 失效依赖清理。原配置、网格、速度表和运行脚本逐字节未变。

删除名单、原行数及逐文件 SHA-256 见 `cleanup_evidence/deletions.json`；静态检查见 `cleanup_evidence/static_preservation.json`。历史发布清单保留在 `provenance/cleanup_origin/`，当前身份见根 `SOURCE_IDENTITY.json`；根 `SHA256SUMS` 对应本次精简交付。

## 当前验证

| 层级 | 本次实际结果 | 结论边界 |
|---|---|---|
| 编译 | 原版、精简版各从独立空 obj/bin 编译，GNU MPI，exit=0 | 未运行 Intel/Slurm |
| 冷启动冒烟 | 高阶隐式、二阶隐式、二阶显式，每版每方法 2 步，8 ranks，全部 exit=0 | 只将两个迭代上限改为 2，物理参数/分块/门不变 |
| 完整重启 | 两种隐式分别从相同原版两步完整态续算 1 步，每版均成功，累计状态号为 3 | 不包含显式续算或跨分区续算 |
| 数值回归 | 三种冷启动和两种续算的 W/h/b 完整保存态及相应输出逐字节一致 | 仅认证这些相同平台短轨迹 |
| 残差/壁面 | 接受态残差、壁面数据、质量伙伴、原门及参考尺度逐字节一致；混合 PLT 的所有非计时字段一致 | 耗时列显式排除，不据此声称提速 |
| 科学验证 | 未新增 | 不改变原高阶残差平台结论，不认证空间阶、物理解精度或长程收敛 |

冷启动完整文件对照数为高阶隐式 30、二阶隐式 22、二阶显式 9；两种续算分别 30、22。每个对照文件的哈希、混合表排除的计时列见 `cleanup_evidence/runtime_parity.json`。首次受沙箱限制的 MPI 端口绑定失败发生在求解器启动前；允许本机通信后完成上述检查，失败日志亦保留在研究目录。

## 使用与复核

继续使用原 `build.sh`、`run.sh`、`resume.sh`、`cases/` 和超算脚本。数值入口和参数卡保持原样；本次没有自动启动正式长程轨迹。与原版同样，修改正式工况或更换平台需要相应资格检查。

研究目录包含 `baseline/`、本目录 `simplified/`、`cleanup_source.patch`、`simplify.py`、`verify_static.py`、`verify_runtime.py` 和原始验证输出。源补丁可用 `patch -p1` 应用于原远端提交的独立副本；它仅修改 src/Makefile，不包含交付元数据。

本次保留仍有调用的旧入口、诊断、环境变量分支以及未确认的冗余。删除的过程包含旧模块公开符号；独立仓外程序若直接链接这些符号，应保留原版或调整调用。本交付保证范围为本仓库现有入口及已验证配置，不承诺旧辅助库 API 完整兼容。
