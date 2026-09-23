#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
run_root="$(mktemp -d)"
trap 'test -n "$run_root" && rm -rf "$run_root"' EXIT

mkdir -p "$run_root/base" "$run_root/reference" \
  "$run_root/dt0040" "$run_root/dt0020" "$run_root/dt0010"
cp "$repo_root/bin/ugks_solver" "$repo_root/tests/cases/implicit_unsteady_uniform/Config.ini" "$run_root/base/"
python3 "$repo_root/tests/helpers/make_cartesian_mesh.py" "$run_root/base/Mesh2d.dat" \
  --nx 16 --ny 4 --surf "$run_root/base/surf.dat"
python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/base/Config.ini" \
  "TimeIntegration_method='Explicit'" Iflag_InitFlow=0 Iflag_WriteRestart=1 \
  Nx_all=16 Ny_all=4 MAX_ITER=1 Nstep_show=1 Nstep_write=1 dt_global=1.25d-3
if ! (cd "$run_root/base" && mpirun -np 1 ./ugks_solver > solver.log 2>&1); then
  cat "$run_root/base/solver.log" >&2
  exit 1
fi

for name in reference dt0040 dt0020 dt0010; do
  cp "$repo_root/bin/ugks_solver" "$repo_root/tests/cases/implicit_unsteady_fwave/Config.ini" \
    "$run_root/base/Mesh2d.dat" "$run_root/base/surf.dat" \
    "$run_root/base/flow2d.dat" "$run_root/base/vdf2d.dat" "$run_root/$name/"
done
python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/reference/Config.ini" \
  "TimeIntegration_method='Explicit'" dt_global=1.25d-3 MAX_ITER=64 Nstep_show=64 Nstep_write=64
python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/dt0040/Config.ini" \
  dt_global=4.0d-2 MAX_ITER=2 Nstep_show=2 Nstep_write=2 Unsteady_InnerTol=5.0d-8 Unsteady_MaxInner=50
python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/dt0020/Config.ini" \
  dt_global=2.0d-2 MAX_ITER=4 Nstep_show=4 Nstep_write=4 Unsteady_InnerTol=5.0d-8 Unsteady_MaxInner=50
python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/dt0010/Config.ini" \
  dt_global=1.0d-2 MAX_ITER=8 Nstep_show=8 Nstep_write=8 Unsteady_InnerTol=5.0d-8 Unsteady_MaxInner=50

for name in reference dt0040 dt0020 dt0010; do
  if ! (cd "$run_root/$name" && mpirun -np 1 ./ugks_solver > solver.log 2>&1); then
    cat "$run_root/$name/solver.log" >&2
    exit 1
  fi
done
python3 "$repo_root/tests/helpers/check_implicit_unsteady_temporal.py" "$run_root"
