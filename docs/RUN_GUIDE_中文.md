# 配置与复现

所有命令从包根目录执行，算法只读取运行目录的 `Config.ini` 和 `Mesh2d.dat`。`run_case.py` 负责复制输入、清理遗留研究环境变量、设置单线程及单 MPI rank，并拒绝覆盖已有输出。

## 例子

| case | 目的 | 是否重启 | 运行量 |
|---|---|---|---|
| `bgk_saved_step` | C4 N16 BGK 保存态复现 | 是 | 1 个 S2CN2 步 |
| `shakhov_saved_step` | N16、Pr=2/3 保存态复现 | 是 | 1 个 S2CN2 步 |
| `steady_saved_step` | N8、累计1300态上的单修正 | 是 | 1 个隐式修正 |
| `steady_finish` | 复现 I3-R1 终止修复 | 是 | 预计9次真实修正至连续五次达标，上限200 |
| `bgk_cold_smoke` | 原始光滑初态接线 | 否 | 2 个物理步 |
| `shakhov_cold_smoke` | 原始 Shakhov 光滑初态接线 | 否 | 2 个物理步 |
| `steady_cold` | 从原初态独立求稳态 | 否 | 上限1500；本次整理不运行长轨迹 |

前四项包含 `restart.bin` 和 `expected_final.bin`。示例依赖全部在包内。保存态是在本机 GNU Fortran/Open MPI 下产生的原生 stream 二进制；跨体系结构可能有端序／种类布局差异。新编译器数值不逐位一致时脚本会报告比较失败，需另做数值差异分析，不能自动当作算法失败或随意放宽门。

```sh
python3 -B scripts/build.py
python3 -B scripts/run_case.py steady_finish --output runs/steady_finish_01
```

构建脚本在 macOS 自动获取 SDK，运行 `make release FC=mpifort`，继承冻结 `-O3 -march=native -funroll-loops`；因此应在目标机器重新构建。并执行小型终止逻辑测试。GNU/Open MPI 是本包实测配置，原 Makefile 的 Intel 分支未在本包验证。MPI 若因沙盒端口权限退出，属于环境问题，应在具有正常本地 MPI 权限的终端复现，不改变科学参数。

## 不要误读的配置

- `TimeIntegration_method='S2CN2'`：CN2 物理时间推进；`S2Steady`：定常隐式修正。
- `dt_global` 同时确定 S2 面积分的物理区间；定常时也不能拿伪时间随意替代它。
- 当前 `tau0=.1` 定义在 `mod_S2D1.f90`，不是通用由 Config 黏度律控制的 variable-τ 路径。
- `Ma=0` 是沿用的初始化配置字段；实际周期初态采用适配器固定 `U=.6,V=.3,T=1`、`ρ=1+.01 sin(x+y)` 的精确单元平均。不能凭这个 Ma 字段解释当前物理流场。
- `Pr=1` / `2/3` 控制 BGK／Shakhov；CK 必须为1。
- 稳态例的 `Implicit_nCFL=1`、`CFL=.1`，初始化得到伪时间 .05；宏观2次、微观1次 sweep。默认不改这些冻结值。
- 当前 GH16×16 与 N8/16/32/64 限制由适配器直接检查，不能仅改 Nu/Nv 绕过。

## 输出与重启

- `native_state_initial.bin` / `native_state_final.bin`：完整 W/h/b/u/v/weight。
- `native_field.csv`：宏观场与实际热流。
- `native_steps.csv`：物理步门值和分项成本。
- `s2i1_true_residual.csv`、`s2i1_updates.csv`：真实稳态残差与接受修正。
- `ConvergenceStatus.txt`：是否满足稳态门；`steady_finish` 预期 `S2Steady T 9`。
- `launch.json`、`stdout.log`、`result.json`：本次复现的真实命令、返回码和字节比较。

原生保存格式（本机 integer32/real64）按顺序为 nx、ny、nvel、time，再 W、h、b、u、v、weight；数组按 Fortran 列主序。重启目前重新计本地迭代和时间字段，旧时间单独写到 `o1_restart_origin.txt`。不要把新输出编号9解释为整个轨迹只用了9步；本例累计为1309。续行比较必须按原初态尺度换算残差，并继承原守恒基准。

## 封存和重复运行

`python3 -B scripts/verify.py` 不执行数值计算。它检查交付文件与源码身份；构建和运行只生成 `bin/obj/runs`，不会改算法。构建日志在 `verification/build.*`，再次构建会更新这些日志，因此新实验应在整个包的副本中进行，或把其视为工作副本，不要重新封存并冒充原记录。
