# v1.2：方腔近壁 2→3→5 重构

近壁第一层复用 Van Leer，第二层采用紧致三阶 WENO，模板充足后恢复 WENO5-AO；角区按 x/y 两方向分别处理。公共面共用通量，保留原 PR 宏观预测、微观单次前后 LU-SGS 扫描及原数值门。包括四壁方腔接线与 GH28 速度表。

版本基于 v1.1 的已验证方腔开发线，源码和 Makefile 与 R23、R25 的 2→3→5 算例版本完全一致；没有合入 R5–R9 复杂隐式原型，也没有额外叠加 v1.1a 清理改动。原有 Couette 配置保留。

## 已有验证

- N16：冷启动536次迭代，五个连续接受态通过原稳态门，最后选定残差最大值约3.27e-9。
- N32：续行至累计2132次迭代，五个连续接受态通过原稳态门，最后选定残差最大值约7.33e-9。
- 这些结果来自已有封存计算，本次发布不重跑轨迹；不代表一般壁面全空间五阶、所有工况准确性或 Intel/超算新版本验证。
- 本次发布检查见 `docs/v1.2_build_check.txt`；历史检查目录中的报告保留其原版本范围。

## 编译与启动

```bash
make release
mkdir -p runs/cavity_N16
cp examples/cavity_graded235_N16/* runs/cavity_N16/
bash run_graded235.sh "$PWD/runs/cavity_N16"
```

示例为 Kn_HS=0.075、N16×16、GH28×28、8 rank（2×4）、nCFL=10，闭腔比例质量校正开启。MPI槽位不足的本机可显式设置 `MPI_OVERSUBSCRIBE=1`。超算须申请8核并加载匹配的MPI编译环境；Intel构建可用 `make release FC=mpiifort`，本版本未重新验证Intel。

`run_graded235.sh` 显式开启 `S2_NEAR_WALL_VL=1` 与 `S2_NEAR_WALL_GRADED235=1`，后者优先；直接使用旧 `run.sh` 不会自动开启2→3→5。不要只改方法标签而遗漏开关。已有 `CORNER_IMPLICIT_DIAG` 等诊断保持关闭；本版本保留这些休眠入口以维持验证源码身份。

流场、残差与续算频率由 Config.ini 控制。运行目录独立，不覆盖旧结果；稳态残差不是物理时间精度。后续重点为更多壁面工况、边界层验证和跨编译器检查。
