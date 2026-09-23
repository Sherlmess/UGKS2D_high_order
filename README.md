# UGKS2D_high_order — v1.0

**二维宏微观高阶重构＋二阶时间推进，附定常隐式加速。**

[中文阅读与运行入口](README_中文.md) · [v1.0 已实现与后续工作](docs/RELEASE_v1.0_中文.md) · [算法指南](docs/ALGORITHM_GUIDE_中文.md) · [复现指南](docs/RUN_GUIDE_中文.md)

The maintained S2 route uses macro/micro WENO5-AO reconstruction, characteristic UGKS face fluxes, a low-cost S2 collision closure with discrete 2×2 heat-flux coupling, and CN2 physical-time integration. It also includes steady implicit acceleration using the same high-order residual, O1–O3 equivalent optimizations, and the corrected steady termination logic.

## Validated scope

Single MPI rank, uniform Cartesian doubly periodic grids, N=8/16/32/64, GH16×16 velocity quadrature, CK=1 and constant tau=0.1. BGK/Shakhov smooth physical-time tests support second-order time convergence. Full steady implicit agreement and observed convergence-cost gain are limited to the archived N8 BGK test.

This release does not certify general fifth-order full spatial discretization, NS/AP behavior, S2 wall boundaries, stretched meshes, multi-rank S2, or unsteady implicit S2.

## Build and reproduce

GNU Fortran, Open MPI, make and Python 3 are required; macOS builds also need Command Line Tools.

```sh
git clone git@github.com:Sherlmess/UGKS2D_high_order.git
cd UGKS2D_high_order
git checkout v1.0
python3 -B scripts/verify.py
python3 -B scripts/build.py
python3 -B scripts/run_case.py bgk_saved_step --output runs/bgk_check_01
python3 -B scripts/run_case.py steady_finish --output runs/steady_finish_01
```

All required example inputs and expected saved states are included. Rebuild on the target machine; no platform-specific solver executables or object files are distributed. Do not interpret same-platform bytewise regression as a cross-platform reproducibility guarantee.

## Provenance and related repositories

- Native foundation: [implicitUGKS v2.0](https://github.com/Sherlmess/implicitUGKS/tree/v2.0).
- One-dimensional S2: [UGKS_high_order_validation, v0.1.0-s2-cn2-1d](https://github.com/Sherlmess/UGKS_high_order_validation/tree/v0.1.0-s2-cn2-1d/deliverables/S2_1D_MacroMicroHighOrder_CN2_20260921).
- Full source identity and historical validation: `provenance/` and `verification/`.
