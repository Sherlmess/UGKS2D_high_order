#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
run_root="$(mktemp -d)"
trap 'test -n "$run_root" && rm -rf "$run_root"' EXIT

mkdir -p "$run_root/explicit" "$run_root/implicit"
for mode in explicit implicit; do
    if find "$run_root/$mode" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
        echo "refusing to reuse non-empty output directory: $run_root/$mode" >&2
        exit 1
    fi
    cp "$repo_root/bin/ugks_solver" "$repo_root/tests/cases/implicit_unsteady_uniform/Config.ini" \
        "$run_root/$mode/"
    python3 "$repo_root/tests/helpers/make_cartesian_mesh.py" "$run_root/$mode/Mesh2d.dat" \
        --nx 8 --ny 6 --surf "$run_root/$mode/surf.dat"
    python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/$mode/Config.ini" \
        Iflag_InitFlow=0 Iflag_WriteSurf=1 Iflag_WriteRestart=1 \
        Nx_all=8 Ny_all=6 MAX_ITER=1 Nstep_show=1 Nstep_write=1 \
        Ma=0.5d0 Re=100.0d0 Tinf=1.0d0 Tw=1.1d0 \
        dt_global=1.0d-6 CFL=0.5d0 \
        "BC_Type=4 6 2 4" FluxWall_type=2 ACC=0.6d0 \
        Quadrature_method=2 Nu=16 Nv=16 NP_X=1 NP_Y=1 \
        Unsteady_Theta=0.5d0 Unsteady_InnerTol=1.0d-10 Unsteady_MaxInner=80
done

python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/explicit/Config.ini" \
    "TimeIntegration_method='Explicit'"
python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/implicit/Config.ini" \
    "TimeIntegration_method='ImplicitUnsteady'"

for mode in explicit implicit; do
    if ! (cd "$run_root/$mode" && mpirun -np 1 ./ugks_solver > solver.log 2>&1); then
        cat "$run_root/$mode/solver.log" >&2
        exit 1
    fi
done

python3 "$repo_root/tests/helpers/check_implicit_unsteady_acc06.py" "$run_root"
