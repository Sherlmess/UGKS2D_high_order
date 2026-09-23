#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
run_dir="$(mktemp -d)"
trap 'test -n "$run_dir" && rm -rf "$run_dir"' EXIT
cp "$repo_root/bin/ugks_solver" "$repo_root/tests/cases/implicit_uniform/Config.ini" "$run_dir/"
python3 "$repo_root/tests/helpers/make_cartesian_mesh.py" "$run_dir/Mesh2d.dat" --nx 4 --ny 4 --surf "$run_dir/surf.dat"
(cd "$run_dir" && mpirun -np 1 ./ugks_solver)
python3 "$repo_root/tests/helpers/check_implicit_uniform.py" "$run_dir"
