#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
run_dir="$(mktemp -d)"
cleanup() {
    find "$run_dir" -mindepth 1 -maxdepth 1 -type f -delete
    rmdir "$run_dir"
}
trap cleanup EXIT
cp "$repo_root/bin/ugks_solver" "$repo_root/tests/cases/implicit_uniform/Config.ini" "$run_dir/"
python3 "$repo_root/tests/helpers/patch_config.py" "$run_dir/Config.ini" \
  "TimeIntegration_method='ImplicitSteady'" Quadrature_method=2 Nu=16 Nv=16 \
  MAX_ITER=5 Implicit_MaxIter=5 Iflag_InitFlow=0
python3 "$repo_root/tests/helpers/make_cartesian_mesh.py" "$run_dir/Mesh2d.dat" \
  --nx 4 --ny 4 --surf "$run_dir/surf.dat"
(cd "$run_dir" && mpirun -np 1 ./ugks_solver | tee solver.log)
python3 "$repo_root/tests/helpers/check_implicit_gh16_uniform.py" "$run_dir"
