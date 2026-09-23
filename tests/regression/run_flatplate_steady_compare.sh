#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
work_root="${RESULT_ROOT:-$(mktemp -d)}"
if [[ -z "${RESULT_ROOT:-}" ]]; then
  trap 'test -n "$work_root" && rm -rf "$work_root"' EXIT
fi
max_iter="${MAX_ITER_LIMIT:-200000}"
show_interval="${SHOW_INTERVAL:-100}"

run_case() {
  local mode="$1" np_x="$2" np_y="$3" name="$4"
  local run_dir="$work_root/$name"
  if [[ -e "$run_dir" ]]; then
    echo "refusing to reuse existing run directory: $run_dir" >&2
    exit 1
  fi
  mkdir -p "$run_dir"
  cp "$repo_root/bin/ugks_solver" "$repo_root/Config.ini" "$repo_root/Mesh2d.dat" "$run_dir/"
  python3 "$repo_root/tests/helpers/patch_config.py" "$run_dir/Config.ini" \
    "TimeIntegration_method='$mode'" Quadrature_method=1 Iflag_WriteRestart=0 \
    CFL=0.5d0 EPS=1.0d-8 Implicit_ResidualTol=1.0d-8 "MAX_ITER=$max_iter" \
    "Implicit_MaxIter=$max_iter" "Nstep_show=$show_interval" "Nstep_write=$show_interval" \
    "NP_X=$np_x" "NP_Y=$np_y"
  (cd "$run_dir" && mpirun -np "$((np_x * np_y))" ./ugks_solver)
  python3 "$repo_root/tests/helpers/check_steady_convergence.py" "$mode" "$run_dir" --tol 1.0e-8
  local flow_file
  flow_file="$(python3 -c 'import pathlib,re,sys; files=list(pathlib.Path(sys.argv[1]).glob("Flow2d_[0-9]*.plt")); print(max(files,key=lambda p:int(re.search(r"_(\d+)\.plt$",p.name).group(1))))' "$run_dir")"
  cp "$flow_file" "$work_root/$name.plt"
  cp "$run_dir/Surf_CpCfCh.plt" "$work_root/$name-surface.plt"
}

run_case Explicit 2 2 explicit4
run_case ImplicitSteady 2 2 implicit4
run_case ImplicitSteady 1 1 implicit1
python3 "$repo_root/tests/helpers/compare_flow_fields.py" "$work_root/explicit4.plt" "$work_root/implicit4.plt" \
  --reference-surface "$work_root/explicit4-surface.plt" --candidate-surface "$work_root/implicit4-surface.plt"
python3 "$repo_root/tests/helpers/compare_flow_fields.py" "$work_root/implicit1.plt" "$work_root/implicit4.plt" \
  --reference-surface "$work_root/implicit1-surface.plt" --candidate-surface "$work_root/implicit4-surface.plt"
echo "PASS flat-plate explicit/implicit and one/four-rank comparisons"
