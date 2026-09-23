#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
run_root="$(mktemp -d)"
trap 'test -n "$run_root" && rm -rf "$run_root"' EXIT
mkdir -p "$run_root/one_rank" "$run_root/two_rank"
for name in one_rank two_rank; do
  cp "$repo_root/bin/ugks_solver" "$repo_root/tests/cases/implicit_unsteady_uniform/Config.ini" "$run_root/$name/"
  python3 "$repo_root/tests/helpers/make_cartesian_mesh.py" "$run_root/$name/Mesh2d.dat" \
    --nx 4 --ny 4 --surf "$run_root/$name/surf.dat"
done
python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/two_rank/Config.ini" NP_X=2 NP_Y=1
if ! (cd "$run_root/one_rank" && mpirun -np 1 ./ugks_solver > solver.log 2>&1); then
  cat "$run_root/one_rank/solver.log" >&2
  exit 1
fi
if ! (cd "$run_root/two_rank" && mpirun -np 2 ./ugks_solver > solver.log 2>&1); then
  cat "$run_root/two_rank/solver.log" >&2
  exit 1
fi
python3 "$repo_root/tests/helpers/check_implicit_unsteady_uniform.py" "$run_root/one_rank"
python3 "$repo_root/tests/helpers/check_implicit_unsteady_uniform.py" "$run_root/two_rank"
