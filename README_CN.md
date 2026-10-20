# 当前版本 IO-R3：速度表来源与超算启动修复

先读 [VELOCITY_SOURCE_FIX_CN.md](VELOCITY_SOURCE_FIX_CN.md)。本包以配置的公共外部速度表作为高低阶共同来源，使用srun/pmi2，必须在新目录编译。下文IO-R2/IO-R1为继承说明。

# 当前版本 IO-R2：残差增加 c_minus_1

本次新增列、输出时刻和验证见 [RESIDUAL_COLUMNS_CN.md](RESIDUAL_COLUMNS_CN.md)。本页下方IO-R1说明作为继承配置说明；算法与输出频率均未改变。

# Ma3 Couette：三方法超算包 IO-R1（关闭逐步审计快照）

本包基于静止壁修复版本，新增仅影响写文件的输出控制。PR、LU-SGS、通量/碰撞、重构、接受/拒绝门及每步真实残差/连续五态收敛计算保持不变。旧超算包和正在运行的作业未改。

必须在全新目录重新编译本包，旧二进制不能解析新Config参数。源代码测试见 OUTPUT_TESTS_CN.md，具体频率见下表。

## 1. 共同工况

- Nx=8、Ny=64；NP_X=1、NP_Y=8，共8个MPI进程，单节点、每进程1线程。
- 同一 Mesh2d.dat、同一半区速度表 velocity_half32_scale2.dat，32×32=1024个速度节点。
- 初始气体静止；上下壁 x 速度为 -3.872983346207417、+3.872983346207417。Ma=0 表示初始气体，不是取消Ma3移动壁。
- Tinf=Tw=2，原生温标 p=rho*T/2；theta=1。ck=1、Pr=2/3、omega=.5、mu_ref=.05637282460904075。
- x周期，y等温漫反射，FluxWall_type=1、CouetteUseW2BWall=true、UnifiedPhysicalOperators=true。
- 相同物理CFL=.1，TimeStep_method=1。dt_global=.00025是既有启动占位，实际dt按CFL重算，不是固定.00025。
- 隐式nCFL=1000、PR2、微观一前一后单扫、alpha原实现、零重试。仅封闭Couette隐式启用原R6比例质量校正；显式关闭。没有新增约束。

## 2. 完整参数卡与停止条件

|案例|完整配置|方法|最大步数/迭代|输出|
|---|---|---|---:|---|
|高阶隐式|cases/high_implicit/Config.ini|HighOrderImplicit|3000|屏幕100，流场100|
|二阶隐式|cases/low_implicit/Config.ini|ImplicitSteady|3000|屏幕100，流场100|
|二阶显式|cases/low_explicit/Config.ini|Explicit|200000|屏幕100，流场1000|

两隐式保留原残差1e-8、宏微1e-12合同及连续五个真实接受态通过后停止。二阶显式保留其原生更新残差EPS=1e-12后停止。显式更新残差与隐式物理残差不是同一量，因此先对比充分收敛后的速度、温度、热流、壁剪切等，不能直接按结束墙钟宣称同精度加速。达到迭代上限、时间上限或exit=0都不自动证明稳态。

本次目标为定常比较，没有新设t=10硬截止。200000是显式步数上限，不承诺相应物理终时；物理dt会随状态变化。若仍未收敛，应使用完整保存态续算，不能拿未收敛末态判断方法准确性。root Config.ini是高阶隐式参数卡的副本；提交脚本只读取cases目录。examples仍保留旧包的三步演示，不用于此次正式提交。

## 3. 上传、编译和提交

解压至新目录，例如 /work/home/sjzhang/Linear_UGKS/Couette_3Methods_IO_R1，不覆盖旧作业目录。在新包根目录执行：

```sh
bash build_hpc.sh
bash submit_high_implicit.sh
bash submit_low_implicit.sh
bash submit_low_explicit.sh
```

三个命令各提交一个8核作业；也可以只提交其中一个。无需chmod，直接bash。不要用sbatch提交submit_*.sh：它们是提交前冻结输入的包装器。

build_hpc.sh复用已有站点模块 compiler/intel/2021.3.0 与 mpi/intelmpi/2021.3.0，使用mpiifort和icc。job.slurm保留wzhcnormal分区及 /opt/gridview/slurm/lib/libpmi2.so，申请1节点8核24小时。若站点模块或队列限制有变化，请按站点实际信息调整；本机未验证Intel/调度器。也可加载环境后手动：

```sh
make -j4 FC=mpiifort COMPILER=intel CC=icc bin/ugks_solver
```

每次提交会创建不重复的 runs/<方法>_<时间>_<随机串>/，先复制Config.ini、网格、速度表、当前可执行文件和job.slurm，再从该目录sbatch。之后改cases配置或重新编译新包，都不会改变已提交作业的副本；不要手动修改runs里的输入。INPUT_SHA256SUMS在作业开始前检查冻结身份。不要在旧代码目录make clean，也不要替换旧作业中的可执行文件。

## 4. 运行时间控制和存储

隐式R6原默认230秒软停止已通过环境变量R5_WALL_LIMIT=82800（23小时）覆盖；数值公式不变，留1小时给收尾。触发软停止仅表示受控保存，不是收敛。显式不使用这个R6控制，由其原生EPS、MAX_ITER和Slurm24小时上限控制。

隐式可以在对应run目录创建STOP_REQUEST请求受控停止。本包不声称此选项对普通Explicit有效。Slurm强杀可能赶不上最终输出，应保留周期结果和完整续算文件。

### Config.ini 输出频率

三份完整Config均已有这些键：

```fortran
Nstep_show    = 100    ! 屏幕及原生残差PLT间隔；必须>=1
Nstep_write   = 100    ! 隐式正常流场间隔；显式配置为1000
Nstep_restart = 100    ! 覆盖最新完整续算；0=只在正常结束/受控停止保存
Nstep_monitor = 1      ! 接受态标量残差、质量/校正/壁面CSV；必须>=1
Nstep_audit   = 0      ! 历史state/restart/profile快照；0=关闭
```

默认不生成逐步 `state_*`、`restart_*`、`profile_*`、`reference_<step>.in`。需要审计时可设Nstep_audit为100等正整数；历史完整快照在间隔及初末/强制事件写出。不要为正式长算设成1。

保留正常流场Flow2d*.plt、Surf.plt/surf.dat与残差文件。必要完整续算每100步覆盖同一套flow2d.dat、vdf2d.dat、closed_mass_rank*.bin和residual_reference.in；显式不需要后两类。最终正常退出/隐式受控停止时，不论是否到间隔都保存最终流场和最新完整续算。Iflag_WriteRestart继续保持1。

每步仍计算真实残差并判定收敛，Nstep_monitor只控制序列化；默认1保留完整标量历史，数据量远小于微观数组。Nstep_show原生残差可能对应更新前状态，隐式最新接受态看sparse_residual.csv。Nstep_monitor>1时cost_steps只保留抽样步的分项耗时，不能把抽样行之和当成完整运行成本。最终残差行强制保留。

8×64 / 32×32速度配置下，当前隐式最新完整续算一套实测约26.05 MiB，大小不随迭代数增长。四步完整目录约27 MiB。按默认流场/标量频率，3000步预计为几十到约100 MiB/条（是文件尺寸外推，不是长程运行实测），不再是77 GiB；实际以输出文件统计为准。低阶显式200000步的历史正常流场也会按1000步间隔累积，通常仍为百MiB量级。

沿用原生续算格式与覆盖写方式，不新增断电原子事务。写文件过程中强杀/掉盘不保证最新一套完整；优先用受控停止并待写完。不要在旧作业目录替换源码、可执行文件、Config或删除旧文件。

## 5. 需要返回的结果

优先返回每个run的stdout.log、slurm日志、job_status.txt、Config.ini、INPUT_SHA256SUMS、Flow2dC.plt/最终Flow2d_*.plt、Surf.plt/surf.dat、ResidualImplicit.plt或原生显式残差文件。隐式另返回R5_status.txt、sparse_residual.csv、consecutive_states.csv、light_steps.csv、cost_steps.csv、wall_history.csv。

完整续算保留flow2d.dat、vdf2d.dat、全部closed_mass_rank*.bin（隐式）、residual_reference.in（隐式）、配置、网格及速度表。程序在每次完整续算写出时同步保存residual_reference.in，不再依赖历史reference快照。job脚本仍兼容审计模式的旧reference复制；强杀不保证收尾。不要只保留图片或删掉伙伴文件。续算仍需相同分块，操作说明见README_SOURCE_CN.md和原resume.sh；当前三个提交入口是冷启动，不能用于直接续算。

## 6. 本次检查范围

GNU编译通过。真实8rank短测：三方法旧/新稀疏/审计输出的四步W/h/b逐位一致；两隐式完整质量伙伴逐位一致；两隐式连续4步与2+2续算逐位一致；间隔外终态保存、零追加步受控停止保存与负频率初始化前拒绝通过。真实残差与壁面数据在共同输出态逐值一致。未做长程、精度认证、Intel/超算实测或正式作业提交。

原完整Config物理设置、8rank分块、两隐式3000上限与显式200000上限、24小时提交模板和23小时隐式软停止保持。提交方法不变：bash build_hpc.sh，然后bash submit_high_implicit.sh / submit_low_implicit.sh / submit_low_explicit.sh。
