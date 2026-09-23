#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
work_root="${RESULT_ROOT:-$(mktemp -d)}"
if [[ -z "${RESULT_ROOT:-}" ]]; then
  trap 'test -n "$work_root" && rm -rf "$work_root"' EXIT
fi
max_iter="${MAX_ITER_LIMIT:-50000}"
residual_tol="${RESIDUAL_TOL:-1.0e-8}"
mach="${MACH:-0.5d0}"
read -r -a accommodations <<< "${ACCOMMODATIONS:-0.0 0.5 1.0}"

run_case() {
  local mode="$1" accommodation="$2" label="$3"
  local run_dir="$work_root/$label"
  if [[ -e "$run_dir" ]]; then
    echo "refusing to reuse existing run directory: $run_dir" >&2
    exit 1
  fi
  mkdir -p "$run_dir"
  cp "$repo_root/bin/ugks_solver" "$repo_root/tests/cases/implicit_uniform/Config.ini" "$run_dir/"
  python3 "$repo_root/tests/helpers/patch_config.py" "$run_dir/Config.ini" \
    "TimeIntegration_method='$mode'" TimeStep_method=2 Iflag_InitFlow=0 \
    Ma="$mach" Re=100.0d0 AOA=0.0d0 Tinf=1.0d0 Tw=1.0d0 \
    CFL=0.5d0 EPS="$residual_tol" Implicit_ResidualTol="$residual_tol" \
    MAX_ITER="$max_iter" Implicit_MaxIter="$max_iter" Nstep_show=500 Nstep_write="$max_iter" \
    Quadrature_method=2 Nu=16 Nv=16 Iflag_WriteSurf=1 Iflag_WriteRestart=0 \
    Nx_all=8 Ny_all=6 NP_X=1 NP_Y=1 "BC_Type=4 6 2 4" BG_geometry=0 \
    FluxWall_type=2 ACC="$accommodation"
  python3 "$repo_root/tests/helpers/make_cartesian_mesh.py" "$run_dir/Mesh2d.dat" \
    --nx 8 --ny 6 --surf "$run_dir/surf.dat"
  (cd "$run_dir" && mpirun -np 1 ./ugks_solver | tee solver.log)
  python3 "$repo_root/tests/helpers/check_steady_convergence.py" "$mode" "$run_dir" --tol "$residual_tol"
  local flow_file
  flow_file="$(python3 -c 'import pathlib,re,sys; files=list(pathlib.Path(sys.argv[1]).glob("Flow2d_[0-9]*.plt")); print(max(files,key=lambda p:int(re.search(r"_(\d+)\.plt$",p.name).group(1))))' "$run_dir")"
  cp "$flow_file" "$work_root/$label.plt"
  cp "$run_dir/Surf_CpCfCh.plt" "$work_root/$label-surface.plt"
}

for accommodation in "${accommodations[@]}"; do
  suffix="${accommodation//./p}"
  run_case Explicit "$accommodation" "explicit-acc-$suffix"
  run_case ImplicitSteady "$accommodation" "implicit-acc-$suffix"
  python3 "$repo_root/tests/helpers/compare_flow_fields.py" \
    "$work_root/explicit-acc-$suffix.plt" "$work_root/implicit-acc-$suffix.plt" \
    --reference-surface "$work_root/explicit-acc-$suffix-surface.plt" \
    --candidate-surface "$work_root/implicit-acc-$suffix-surface.plt"
  if [[ "$accommodation" == "0" || "$accommodation" == "0.0" ]]; then
    python3 "$repo_root/tests/helpers/check_specular_surface.py" \
      "$work_root/explicit-acc-$suffix-surface.plt"
    python3 "$repo_root/tests/helpers/check_specular_surface.py" \
      "$work_root/implicit-acc-$suffix-surface.plt"
  fi
done
echo "PASS flat-plate Maxwell accommodation explicit/implicit comparisons"
