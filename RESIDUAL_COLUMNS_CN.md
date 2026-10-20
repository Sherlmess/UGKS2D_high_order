# IO-R2：在原残差文件追加质量比例校正系数

本包继承 IO-R1 的三方法配置和低输出策略。仅追加标量输出；不改变PR、微观前后单扫、比例校正、残差计算、收敛门、物理参数或输出频率。不修改旧包及已提交作业。

## 新列与文件

| 方法/文件 | 新列 |
|---|---|
| 高阶/二阶隐式：ResidualImplicit.plt | c_minus_1、mass_rescale_enabled |
| 高阶/二阶隐式：sparse_residual.csv | c_minus_1、mass_rescale_enabled |
| 二阶显式：Residuals.plt 和 Residuals_detail.plt | c_minus_1、mass_rescale_enabled |

两列均追加在最后，旧列位置不变。c_minus_1 是有符号的 c-1，保留约17位有效数字；正值表示这次校正增加质量，负值表示减少质量。mass_rescale_enabled=1 表示启用，0表示关闭。显式两列为0，表示未进行校正，不表示显式物理残差为零。

隐式系数直接读取原算法已经使用的mass_scale，不从校正后质量重新计算，不增加MPI归约、RHS或状态更新。它表示到达当前接受态所应用的最近一次校正，而非多步累计值。

## 输出时刻和阅读方法

- ResidualImplicit.plt：频率仍为Nstep_show；iter为本段局部接受步数。新增c_minus_1对应该次接受更新。原残差列保留原有口径，循环内可能来自更新前状态，末态行可重新评估，不得据此把每行所有量强行视作同一状态的物理残差。
- sparse_residual.csv：频率仍为Nstep_monitor，cumulative_state给出累计接受态编号；这是比较隐式当前接受态原始物理残差和c_minus_1的优先文件。正常末态强制输出规则不变。
- 初始/重启基态尚未执行本段新校正时，c_minus_1=0是初始化占位，不代表上一个运行最后一次校正已知或为零；从本段首个真实接受更新后读取实测系数。
- 二阶显式残差频率仍为Nstep_show；它没有质量比例校正，因此c_minus_1恒为0。
- 若需要每步的原生PLT，设Nstep_show=1。默认隐式Nstep_monitor=1已经保留每步的接受态轻量CSV。Nstep_audit继续为0，不产生逐步大数组审计。
- 绘对数收敛图时可使用abs(c_minus_1)，原文件保留正负号；零值不要直接取log。

没有新运行输出文件，每条已有残差记录只多两个标量。不要在旧残差文件后直接拼接新列格式；在新运行目录使用本包，旧续算输入可以按原说明复制进入。

## 编译与提交

在新目录解压后：

```sh
cd Couette_3Methods_IO_R2
bash build_hpc.sh
bash submit_high_implicit.sh
bash submit_low_implicit.sh
bash submit_low_explicit.sh
```

三个提交脚本与IO-R1逐字节一致，各自冻结输入和可执行文件到独立运行目录。三份完整Config.ini仍位于cases下，参数未变。正在运行的旧作业不会自动增加此列，也不应覆盖其可执行文件。

## 本轮验证

- GNU Fortran/MPI本机编译通过；未新增Intel/超算编译认证。
- 真实8rank：三方法各旧/新4步，共24步；W/h/b完整重启逐字节相同，两隐式质量账本伙伴也相同。
- 两隐式各从已保存2步前缀续行2步，与连续4步逐字节相同。
- 新c_minus_1逐行等于原light_steps.csv实测c_minus1；显式为0且开关列为0。
- 原残差科学字段不变（计时字段不要求一致）；运行目录文件集合不变，无逐步state/restart/profile审计文件。
- 209个旧包文件哈希不变。只改3个输出相关Fortran文件及Makefile依赖，不改数值内核。
- 本轮没有长程收敛、物理准确性或性能优势的新结论。

工程记录：最初日志重定向使用了错误工作目录，被沙箱拒绝，未执行构建；首次源码编译暴露文本替换范围过宽的问题，已修正并最终编译通过。首次MPI启动因沙箱端口限制在求解器运行前退出，保留日志；获得本机通信权限后完成上述短测。无数值拒绝或参数重试。

成功MPI短测进程墙钟合计 14.051 秒。
