#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
run_dir="$(mktemp -d)"
trap 'test -n "$run_dir" && rm -rf "$run_dir"' EXIT
cp "$repo_root/bin/ugks_solver" "$repo_root/Config.ini" "$repo_root/Mesh2d.dat" "$run_dir/"
python3 "$repo_root/tests/helpers/patch_config.py" "$run_dir/Config.ini" \
  Quadrature_method=1 MAX_ITER=3 Nstep_show=1 Nstep_write=3 NP_X=2 NP_Y=2
(cd "$run_dir" && mpirun -np 4 ./ugks_solver)
if [[ "${RECORD_BASELINE:-0}" == "1" ]]; then
  python3 "$repo_root/tests/helpers/check_explicit_baseline.py" "$run_dir" \
    "$repo_root/tests/regression/explicit_structured_3step.json" --record
else
  python3 "$repo_root/tests/helpers/check_explicit_baseline.py" "$run_dir" \
    "$repo_root/tests/regression/explicit_structured_3step.json"
fi
