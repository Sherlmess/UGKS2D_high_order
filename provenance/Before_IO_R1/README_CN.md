# Ma3 Couette：高阶隐式／二阶隐式／二阶显式超算比较

本包冻结于静止壁修复版本 S2D_Implicit_O13_R6_VLWall_R1。`src/`、Makefile 与该包逐字节相同，本轮仅配置和提交接线；没有新增数值验证或改算法。请上传到全新目录，不覆盖正在运行的代码目录。

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

解压至新目录，例如 /work/home/sjzhang/Linear_UGKS/Couette_3Methods_VLWall_R1，不覆盖旧作业目录。在新包根目录执行：

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

当前隐式每个接受态保存state/restart/profile审计文件（源码继承）：按本机同网格文件实测合计约26.2MiB/步，3000步约76.75GiB/条；两个隐式上限约153.5GiB，另加流场等。Nstep_write不控制这些继承审计文件。本轮没有更改输出算法或自动清理，确保磁盘额度足够。正式源码及历史结果保持原样。

## 5. 需要返回的结果

优先返回每个run的stdout.log、slurm日志、job_status.txt、Config.ini、INPUT_SHA256SUMS、Flow2dC.plt/最终Flow2d_*.plt、Surf.plt/surf.dat、ResidualImplicit.plt或原生显式残差文件。隐式另返回R5_status.txt、sparse_residual.csv、consecutive_states.csv、light_steps.csv、cost_steps.csv、wall_history.csv。

完整续算保留flow2d.dat、vdf2d.dat、全部closed_mass_rank*.bin（隐式）、residual_reference.in（隐式）、配置、网格及速度表。job.slurm在正常/失败退出后尽可能复制最后reference到residual_reference.in；强杀不保证执行到该处。不要只保留图片或删掉伙伴文件。续算仍需相同分块，操作说明见README_SOURCE_CN.md和原resume.sh；当前三个提交入口是冷启动，不能用于直接续算。

## 6. 本次检查范围

三份参数卡逐字段比较：除方法、比例校正、步数上限和输出频率外一致；网格/速度表逐字节相同。shell语法和隔离假sbatch/mpirun检查通过，验证了8rank启动、冻结输入、退出码及残差参考传递。假运行中macOS的stack限制由测试桩代替，真实脚本仍执行ulimit -s unlimited。

数值源码未变，沿用上一包必要短测证据；本轮没有编译或执行求解器、没有提交真实作业。高阶与低阶离散平衡构造本身仍各用自己的实现，物理输入相同不表示完整离散算子相同。本包用于你实际计算后验证，不预先认定三条长程解通过。
