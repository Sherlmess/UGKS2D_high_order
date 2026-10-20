# IO-R3：高阶入口统一使用配置选择的外部速度表

## 已定位的问题

上传 high_implicit_20261019T154142_0MOUes 作业46910769已经成功进入Fortran主程序、读取参数和8rank分块。失败为主动MPI_Abort(408)：native velocity table differs from frozen R1 table。尚未开始接受步。不是此前Hydra段错误，也不是c_minus_1输出错误。

配置为StructuredVelocityFile='velocity_half32_scale2.dat'、half32/scale2、Nx8/Ny64、NP_X1/NP_Y8。上传速度表与IO-R2交付表逐字节相同；参数卡也未改变物理配置。

原生初始化先读入外部节点/权重，并由rank0广播。高阶入口又调用init_velocity，按内置原始表使用sqrt/exp重新计算有效权重，要求新算数组与外部文件完全相等。这混用了两种权威来源，编译器数学库、常量折叠及运算顺序差异可触发拒绝。上传日志未保存逐点差值，不能给出Intel实际ULP差，也不能断言仅某一个exp值造成失败。

## 最小修改

仅修改src/mod_S2NativeAdapter.f90：

- 显式配置外部表时，高阶入口以原生已经读入/验证/广播的uSpace、vSpace、weight为来源，检查uVel、vVel、wVel与之完全一致。
- 未配置外部表时，保留原init_velocity生成表及完全一致检查。
- 保留原文件维度、有限性、正权重、张量布局、中心对称检查；入口另显式检查有限性与正权重。
- 不修改任何节点或权重，不覆盖读取数组，不改变CFL、碰撞、重构、隐式扫描、质量校正、残差及科学门。
- 这是外部表来源合同修复：不再要求显式选择的合法外部表必须与另一份内置表逐位相同。输入文件身份继续由作业INPUT_SHA256SUMS固定。

同时继承启动修复job.slurm=srun --mpi=pmi2；提交目录仍隔离。

## 本机验证与限制

1. GNU Fortran/MPI release编译通过。
2. 8rank三方法各旧/新4步，完整W/h/b及隐式质量伙伴逐字节一致；原接受态残差、c-1台账逐值一致。
3. 局部输入合同复现：只将外部文件每个权重上移一个浮点单位并保持对称，旧版在408停止，修复版完成一步。此为诊断数据，不进入正式参数卡。
4. 负权重仍在原读取入口205拒绝；关闭外部文件后的内置生成路径完成一步。
5. 212个旧包文件起止哈希不变；正式配置/网格/速度表与IO-R2一致。

没有Intel编译器和Slurm本地实测。本包支持修复后的本机短步一致性，不能宣布本次超算或长程物理计算已通过。

## 超算操作

上传并在新目录解压：

```bash
unzip Couette_3Methods_IO_R3_source.zip
cd Couette_3Methods_IO_R3
bash build_hpc.sh
bash submit_smoke.sh high_implicit
```

本次修改了源码，需要重新编译。短测最多两步，15分钟、8rank；使用正式高阶参数卡的独立副本，只改变步数上限。job_status.txt中exit=0并且stdout.log确认实际完成两步后，再提交：

```bash
bash submit_high_implicit.sh
# 其他两种方法如需使用统一新包，也可提交：
bash submit_low_implicit.sh
bash submit_low_explicit.sh
```

不要重启旧失败run中的旧二进制。此包不会修改旧目录、正在运行作业或原始失败证据。

继承IO-R2的残差c_minus_1、mass_rescale_enabled、低输出频率和完整续算格式；详细解释见RESIDUAL_COLUMNS_CN.md。首次正式配置仍为定常迭代上限，而非t=10硬截止。

本轮成功与预期拒绝的MPI调用墙钟合计14.436秒。
