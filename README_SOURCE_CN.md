# S2D-VLWall-R1：二阶静止壁平衡修复，保留 O13/S2 与 ZHU-R6

本包在 HI1-R1 上修复二阶统一物理算子的静止壁平衡。仅修改 `src/mod_Evolution.f90` 和 Makefile 依赖；高阶生产源码不变。二阶墙面与内部面使用二阶自身已有的离散平衡矩准备，宏观面通量取同一实际微观通量的离散矩。未采用 S2 非线性平衡求解器，未新增状态/修正量投影或质量约束。PR2、微观一前一后、nCFL=1000、alpha=1 保持不变。

适用本包 `UnifiedPhysicalOperators=true`、全局物理 CFL、W2B/type1 平直壁入口。type0/type2、局部物理步长及开放边界层未新增资格。原关闭分支源码保留，但本包不声明其数值回归通过（旧速度表适用范围限制使关闭探针在初始化前拒绝）。现有高阶路径 Ma3 保存结果逐位回归通过。

旧 HI1 文档仅作为来源放在 `provenance/HI1_R1_documents`，当前测试与结论以本目录 TEST_RESULTS_CN.md 为准。

## 方法与相同物理输入

|目录|TimeIntegration_method|推进|比例质量校正|
|---|---|---|---|
|examples/low_explicit|Explicit|原生 Van Leer + CN 碰撞|关闭|
|examples/low_implicit|ImplicitSteady|原生真实残差 + R6 PR(2)/微观单扫|封闭 Couette 显式开启|
|examples/high_implicit|HighOrderImplicit|S2 真实残差 + R6 低阶近似修正|封闭 Couette 显式开启|
|consistency/high_explicit|S2CN2|保留 S2 高阶显式 CN2|关闭|

三份正式 Config.ini 仅方法名称和必要的质量选项不同；隐式参数在显式入口中不参与推进。统一：Ma3 移动壁、初态气体静止、8×64 网格、32×32 同一外部速度表、Pr=2/3、ck=1、相同黏度、同一上下壁温度和速度、全局物理 CFL=0.1。统一开关 UnifiedPhysicalOperators=true 使显式/隐式采用同一离散平衡矩和流动矩准备。隐式 nCFL=1000、MacroSweeps=2、MicroSweeps=1、MaxRetries=0。时间步方法统一为 1（全局 CFL），S2 每次使用原生更新的物理 dt_global；伪时间步单独计算。

本包所有示例上限为 3 步，是易读启动配置，不是正式收敛上限。正式计算前由用户编辑 MAX_ITER / Implicit_MaxIter；本次没有进行正式计算。源代码缺省质量选项为 false；配置仅封闭周期-x/漫反射-y Couette 隐式显式设为 true。开放边界和低阶显式开启此选项会被拒绝。不得用比例校正代替开放边界的物理质量平衡。

## 编译与本机运行

需 Fortran/MPI、C99 编译器和 make；编译/运行不依赖 Python。GNU 本机已验证：

```sh
bash build.sh
cp -R examples/high_implicit run_hi
MPI_OVERSUBSCRIBE=1 bash run.sh run_hi
```

低阶入口将上面目录换为 low_explicit 或 low_implicit。run.sh 从 Config.ini 读取 NP_X×NP_Y，保留原生工作目录和输出文件。直接运行也可以：进入准备好的算例目录，使用绝对路径 `mpirun -np 8 /path/to/bin/ugks_solver`。O13 优化开关在 run.sh 中启用；直接运行时可同样设置 S2_OWNER_RECON、S2_OWNER_COMPACT、S2_FIRST_LAYER_SCOPE、S2_VELOCITY_TILE、S2_WALL_TILE=1；缺省关闭这些优化不会改变方法公式。

调试构建可用 `make debug`。本次额外用 GNU `-O0 -g -fcheck=all -fbacktrace -finit-real=snan -finit-integer=-999999 -ffpe-trap=invalid,zero,overflow` 在独立干净目录检查数组边界及未初始化读。数组临时量警告属于既有重构布局，不是越界或未初始化错误。

## 提交模板与 Intel 边界

submit.slurm 沿用已有 NativeCouette_original_direction/submit.slurm 的真实站点设置：wzhcnormal、8 rank、compiler/intel/2021.3.0 与 mpi/intelmpi/2021.3.0、原 PMI 路径；未添加账号。原模板逐字节保留在 provenance/original_submit.slurm。只有算例选择和新工作目录拷贝部分改为本包入口。

Intel Makefile 分支保留，可在已有环境加载后用 `make FC=mpiifort COMPILER=intel CC=icc` 构建。本机没有验证 Intel、Intel MPI、该超算模块或调度器运行，不保证其可用。用户应按真实站点现状核对模板；本次未提交任何作业。用户准备完编译和参数后，可自行使用 `CASE_NAME=high_implicit sbatch submit.slurm`。续算应从已准备好的新重启目录启动原生 MPI；保留模板站点资源设置即可，不要把冷启动拷贝命令用于续算。

## 完整续算

封闭 Couette 隐式完整续算需同时保留：Config.ini、Mesh2d.dat、velocity_half32_scale2.dat、flow2d.dat、vdf2d.dat、全部 closed_mass_rankN.bin、residual_reference.in。这些伙伴保留 M0、完整 W/h/b halo、累计数值账本和残差参考尺度，不能丢弃后用当前质量或新残差重新定义它们。分区必须保持原 NP_X / NP_Y；本包没有扩展 R6 伙伴文件的重分区能力。

本包交付源码；首次使用请编译。GNU本机验证环境为macOS ARM64；其他平台也必须在实际MPI环境编译。

正常 run.sh 结束后将最后 reference_*.in 复制为 residual_reference.in。例如从运行末态再追加两步：

```sh
bash resume.sh run_hi run_hi_next 2
```

resume.sh 拒绝覆盖旧目录，仅改变初态方式为 2 和追加迭代上限。不要设置 R6_INITIAL_STREAK；它是旧审查工具的受控选项，不是日常续算参数。只读复核可使用 `R6_MODE=preflight bash run.sh prepared_restart`，实际零更新。缺少或不匹配的质量伙伴/参考文件会被程序拒绝。格式与缓存生命周期未变；旧包物理态可读，后续采用修复后的墙面算子。完整续算仍受原精确身份检查约束，跨O0/O2编译有已观察的求积身份拒绝，不能承诺跨编译兼容；详见ROOT_CAUSE_CN.md。

R6 的 230 秒缺省软控制保留。封闭 Couette 隐式可在新工作目录设定适合已获资源的 R5_WALL_LIMIT（秒）；创建 STOP_REQUEST 后，下一更新前受控停止并保存最后接受态。新续算目录只拷贝必要文件，不拷贝旧 STOP_REQUEST。

## 原生输出

- Flow2dC.plt / Flow2d_*.plt：原生 Tecplot 流场，rho/U/V/T_native/p/应力/热流与原有残差字段；p=rho*T_native/2。
- flow2d.dat / vdf2d.dat：原生宏观与微观重启格式；grid2d.grd 保留原生网格输出。高阶隐式使用 R6 完整重启伙伴；高阶显式原有原子重启元数据入口保留。
- Residuals.plt / Residuals_detail.plt：显式更新残差；ResidualImplicit.plt / ConvergenceStatus.txt：原生隐式更新前缓存口径，不应冒充最新接受态的真实残差。
- 封闭 Couette 的 light_steps.csv、sparse_residual.csv、consecutive_states.csv、profiles.csv、wall_history.csv、R5_status.txt 和 reference_*.in 保留 R6 观察与状态编号。light_steps 的 native_residual_state 为更新前状态，最新接受态看 sparse_residual 和累计编号。
- Surf.plt / surf.dat：原生底壁格式；高阶显式另保留顶壁对应文件，R6 wall_history 记录两侧壁量。高阶壁量由 S2 实际面缓存发布；本例封闭隐式最终残差复算后为末态壁量，普通显式表面可能仍是更新前面缓存。不要把它写成每次都已更新的末态壁量。
- 比例校正账本 dM/dPx/dPy/dE 是数值校正增量，不是壁面物理源。

Ma=0 的共同入口跳过旧外流专用 uvw2D.dat 的 Ma 归一化分支，避免除零；正常原生流场、残差与重启格式没有另造。

## 简单检查与范围

见 TEST_RESULTS_CN.md、ROOT_CAUSE_CN.md 和 SOURCE_PATCH.diff。独立解包构建及静止态启动在本机 GNU/MPI 验证；未验证 Intel/超算。用户后续自行进行长程 Couette 对比，本包没有运行新的稳态或精度矩阵。

静止测试关闭了 ClosedCouetteMassRescale，证明修复不依赖比例校正。Ma3 隐式示例继续使用已经授权的 R6 比例校正；它与本次静止平衡修复是两件事。没有声称三步 Ma3 已稳态，也没有消除历史有限扫描的 rank 依赖。

本次改变低阶物理残差，新的低阶长期结果可能不同于旧包。建议低阶正式对照从原初态开始。新包内部连续3步与1+2完整续算逐位一致；旧低阶态若续算，需注明中途更换离散算子，不能当作单一旧算法轨迹。

## 哪些文件只用于测试

`checks/`、`SOURCE_PATCH.diff`、`TEST_RESULTS_CN.md`、`ROOT_CAUSE_CN.md`、`provenance/` 都不参与编译/推进，结果确认后可从运行副本删除。建议保留一份原始源码包和校验和。`examples/` 和 `consistency/` 是参数样例，准备好正式运行目录后可从计算副本移除。

保留全部 `src/`（包括看似检查用的模块）：当前 Makefile 和主程序有编译依赖，不建议直接删除单个 Fortran 文件。`NATIVE_WIRING_CHECKS` / `NATIVE_OPERATOR_ONLY` 是默认关闭的审计入口，正常计算无需开启。`light_steps`/残差参考/closed_mass_rank 文件属于已有隐式观察和完整续算合同，不应当成可随意删去的静止测试产物。
