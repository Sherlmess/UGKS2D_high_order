# 算法接线与源码差异

入口为同一 main_UGKS2D → Solver → Evolution/ComputeResidual → 原生更新分派，不引入新运行框架。高阶显式 S2CN2 仍可单独选择。新 HighOrderImplicit 复用 S2 的初始化、Cartesian halo、面通量与体求积目标，但更新走 R6，不走旧试验性 S2Steady 的固定 tau=0.1 与全局约束方向。本包不将该旧试验入口列为交付选项。

令 F_hi 为 S2 的真实高阶微观流动残差率，C_hi 为其真实体碰撞率。C_hi 使用 advance_step 原来使用的 mh/mb、kx/ky、heat(h−mh,b−mb) 和 s2_tau(W)，既包含体求积，也保持真实物理黏度；墙面来自 S2 multiscale wall；HI1-R1使其中的墙面目标及反射h/b共用体平衡已有的discrete_equilibrium_q_2d，半空间密度/通量矩闭合保持。宏观右端是同一实际面通量的矩。没有用低阶通量或点平衡替换这些真实残差。

宏观 PR(2) 解原低阶近似修正方程，由高阶宏观右端驱动。取得 W_tilde 后，按 R6 用旧态冻结热流构造低阶点平衡 G_lo(W_tilde;q_old)。微观一次前扫、一次后扫使用：

```
R_hi = F_hi + C_hi
C_lo(W; f_old,q_old) = (G_lo(W;q_old) - f_old) / tau_lo(W)
rhs = R_hi + C_lo(W_tilde;f_old,q_old) - C_lo(W_old;f_old,q_old)
D = 1/dt_pseudo + 1/tau_lo(W_tilde) + 原低阶迎风对角
```

代码中 stream_h/b 只存 F_hi，physical_collision_h/b 存 C_hi，low_collision_old_h/b 存被扣除的旧态低阶碰撞。预测低阶碰撞仅贡献增量；在零宏观方向时增量消失，故所求零残差仍由真实 S2 算子决定。碰撞没有重复加入，也没有漏掉。预测态只用于近似修正，不改写接受态的真实 S2 目标。接受态宏观量从更新后的 h/b 恢复；无额外四矩投影、全局 JVP、GMRES、WENO 导数或微观多扫。

物理面通量已经含 dt_physical 和面长，ComputeResidual 后只用 dt_physical 换算成率；dt_pseudo = dt_physical*nCFL/CFL，只出现在修正对角中。统一全局物理 CFL 使 S2 和原生低阶入口使用相同物理时间尺度，没有将 nCFL 写入物理面通量。

低阶统一入口保留 R6 的流动矩投影以及离散 Maxwell/Shakhov 平衡矩投影。低阶显式的旧态与预测态 CN 碰撞也采用同一平衡矩准备，初始化后同步相同的 halo；因此必须比较实际算子，不依赖配置相同这一推断。保存态实测两对算子的面通量、体碰撞和真实残差逐位相同，但显式与隐式推进后的状态不要求相同。

R6 的 moving-wall 近似修正、MPI 交换、原生重启伙伴、原尺度继承和封闭 Couette 比例校正保留。比例校正同时缩放 W/h/b，恢复宏微一致的接受态，并保留 M0 和账本。它不改变开放边界物理算子，开放边界禁止启用。各入口遵守高精度 PI 和 p=rho*T_native/2。

## 最小接线变化

|文件组|作用|
|---|---|
|Globaldata / ImportPara|新方法选择、共同物理准备开关、R6 参数、质量用途限制|
|S2D1 / S2CartesianMPI / S2NativeAdapter|新增方法别名分派；直接读取原 S2 真实碰撞缓存；每个新接受态重新构造原体目标|
|Residual / ImplicitState|共享真实碰撞缓存及原生残差率；高阶只检查矩一致性，低阶沿用 R6 准备|
|ImplicitMacro / Boundary|R6 PR 与壁面映射，记录实际 PR 次数|
|ImplicitMicro / Diagnostics / SteadyImplicit|R6 单扫描、预测碰撞增量减旧碰撞、真实高阶残差诊断；均匀 S2 初态只核查矩，不再次改变其平衡分布|
|InitializeFlow / Interpolation|相同速度表和初始准备；只初始化完整伙伴会保存的未使用墙外状态；跳过随后被壁面单侧公式覆盖的无效 stencil 读|
|Evolution|原低阶真实面通量加 R6 壁面修复；高阶真实计算仍调用原 S2|
|Solve / Output / R5Observe|同一原生循环、R6 生命周期/完整重启、真实壁缓存输出及零更新输出准备；Ma=0 输出防除零|
|NativeWiringChecks|默认关闭的观察工具，记录真实数组和实际扫描；不增加数值修正|

provenance/O13_to_integrated.patch 与 R6_to_integrated.patch 提供逐源码差异；source_identity.json 保留来源 SHA256。HI1-R1的静止壁离散平衡修复仅复用现有四矩Q，详见ROOT_CAUSE_CN.md和HI1_R1_minimal.patch；未改PR/LU-SGS、重构、气体模型或门。MPI 修正方向仍是原 R6 的有限分块扫描，故一步状态依赖分区；同一保存态的真实目标算子保持一致。
