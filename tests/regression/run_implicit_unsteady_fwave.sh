#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
run_root="$(mktemp -d)"
trap 'test -n "$run_root" && rm -rf "$run_root"' EXIT

mkdir -p "$run_root/base" "$run_root/explicit" "$run_root/implicit"
cp "$repo_root/bin/ugks_solver" "$repo_root/tests/cases/implicit_unsteady_uniform/Config.ini" "$run_root/base/"
python3 "$repo_root/tests/helpers/make_cartesian_mesh.py" "$run_root/base/Mesh2d.dat" \
  --nx 16 --ny 4 --surf "$run_root/base/surf.dat"
python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/base/Config.ini" \
  "TimeIntegration_method='Explicit'" Iflag_InitFlow=0 Iflag_WriteRestart=1 \
  Nx_all=16 Ny_all=4 MAX_ITER=1 Nstep_show=1 Nstep_write=1 dt_global=5.0d-3
if ! (cd "$run_root/base" && mpirun -np 1 ./ugks_solver > solver.log 2>&1); then
  cat "$run_root/base/solver.log" >&2
  exit 1
fi

for mode in explicit implicit; do
  cp "$repo_root/bin/ugks_solver" "$repo_root/tests/cases/implicit_unsteady_fwave/Config.ini" \
    "$run_root/base/Mesh2d.dat" "$run_root/base/surf.dat" \
    "$run_root/base/flow2d.dat" "$run_root/base/vdf2d.dat" "$run_root/$mode/"
done
python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/explicit/Config.ini" \
  "TimeIntegration_method='Explicit'"
if ! (cd "$run_root/explicit" && mpirun -np 1 ./ugks_solver > solver.log 2>&1); then
  cat "$run_root/explicit/solver.log" >&2
  exit 1
fi
if ! (cd "$run_root/implicit" && mpirun -np 1 ./ugks_solver > solver.log 2>&1); then
  cat "$run_root/implicit/solver.log" >&2
  exit 1
fi
python3 "$repo_root/tests/helpers/check_implicit_unsteady_fwave.py" "$run_root"
