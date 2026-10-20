# UGKS2D high order — 二维 Couette 高阶隐式阶段版本

**已成功实现高阶隐式二维 Couette 计算，主要流场与同工况二阶隐式、二阶显式结果接近。**

主方法是宏观／微观高阶重构的 S2，复用原 PR 宏观预测及微观单次前后 LU-SGS 扫描。此处为定常隐式求解；不能把伪时间迭代称为物理时间高阶。

## 已实现

- 高阶隐式 `HighOrderImplicit`、二阶隐式 `ImplicitSteady`、二阶显式 `Explicit` 三个入口及同工况配置。
- Ma3 Couette，8×64 网格，32×32 离散速度，1×8 MPI 分块。
- Shakhov／变 tau、离散平衡一致壁面、已有重构性能优化和 MPI 接线。
- 封闭 Couette 隐式按已有比例校正保持质量，残差输出 `c_minus_1`；显式不开启该校正。
- 正常流场、标量残差及必要完整续算，输出频率在 Config.ini；默认关闭逐步大数组审计。
- 公共外部速度表、Intel/C99 构建及 srun/pmi2 超算脚本。

## 本次 Couette 对照结论

已有三条超算结果的输出流场对照：高阶相对二阶显式，全域最大差按二阶场峰值归一化，密度约0.252%、切向速度0.105%、温度0.133%、压力0.099%、qx/qy约0.884%/0.666%。这些是方法间差异，不是相对精确解的误差。

二阶隐式与二阶显式已输出的宏观量／应力／热流最大绝对差约4.4e-11。高阶的法向小速度分量相对差较大，不能把“主要流场接近”解释为所有分量逐位一致。

二阶隐式满足原连续五态残差门；高阶在3000次上限形成残差平台，尚未满足原严格收敛门。二阶显式到200000步上限，未满足自身停止门。残差平台不自动等于流场错误，但也不据流场接近宣称严格稳态或更高物理准确性。详见 [阶段说明](docs/COUETTE_MILESTONE_CN.md) 与 [已有结果分析](docs/couette_results/REPORT_CN.md)。

## 编译与运行

本机 GNU/MPI：

```sh
bash build.sh
```

超算（在新目录中）：

```sh
bash build_hpc.sh
bash submit_high_implicit.sh
bash submit_low_implicit.sh
bash submit_low_explicit.sh
```

各脚本在独立 runs 目录冻结输入。提交模板按现有站点配置提供，应先核对分区与模块。三个完整配置在 `cases/`；根 Config.ini 是高阶隐式副本。详细温标、参数、输出与续算说明见 [运行指南](README_CN.md)。

## 发布范围

源码来自已交付 IO-R3；本次整理不改变数值算法。清单见 SOURCE_IDENTITY.json，重新构建／最小启动结果见 RELEASE_VERIFICATION.json。旧版本仍可通过 Git 历史和 v1.0 标签访问。

`checks/`、`consistency/` 是定向检查输入；`examples/` 是短测配置；`provenance/` 与 docs 是来源／证据；正式求解需 src、Makefile、配置、网格、速度表和相应运行脚本。src 中诊断模块仍在依赖图中，不应随意删除。

## 后续

- 完成四壁、角点及左右单侧高阶模板，接入 Lim Kn=0.075 方腔；当前 S2 壁面入口仍是 x 周期／y 两壁，不能只改四壁配置就运行。
- 继续区分流场稳定与严格残差达标；评价高阶收益需独立参考或网格比较。
- 新平台／新工况需必要短测；本发布不新增空间阶、NS/AP或普适收敛认证。
