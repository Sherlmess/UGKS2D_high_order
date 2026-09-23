#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
result_root="${RESULT_ROOT:-$(mktemp -d)}"
mkdir -p "$result_root"
csv="$result_root/ncfl_sweep.csv"
max_iter="${MAX_ITER_LIMIT:-200000}"
show_interval="${SHOW_INTERVAL:-100}"
printf '%s\n' 'method,configured_nCFL,min_effective_nCFL,total_retries,converged,iterations,evolution_calls,pr_seconds,lusgs_seconds,correction_exchange_seconds,total_seconds,peak_memory_mb,speedup' > "$csv"

run_one() {
  local method="$1" ncfl="$2" name="$3"
  local run_dir="$result_root/$name"
  if [[ -e "$run_dir" ]]; then
    echo "refusing to reuse existing run directory: $run_dir" >&2
    exit 1
  fi
  mkdir -p "$run_dir"
  cp "$repo_root/bin/ugks_solver" "$repo_root/Config.ini" "$repo_root/Mesh2d.dat" "$run_dir/"
  python3 "$repo_root/tests/helpers/patch_config.py" "$run_dir/Config.ini" \
    "TimeIntegration_method='$method'" "Implicit_nCFL=$ncfl" Quadrature_method=1 \
    CFL=0.5d0 Iflag_WriteRestart=0 "Nstep_show=$show_interval" "Nstep_write=$max_iter" NP_X=2 NP_Y=2 \
    EPS=1.0d-8 Implicit_ResidualTol=1.0d-8 "MAX_ITER=$max_iter" "Implicit_MaxIter=$max_iter"
  (cd "$run_dir" && /usr/bin/time -l mpirun -np 4 ./ugks_solver > solver.log 2> time.log)
  python3 "$repo_root/tests/helpers/check_steady_convergence.py" "$method" "$run_dir" --tol 1.0e-8
}

run_one Explicit 0 explicit
explicit_total="$(awk '/ real /{print $1; exit}' "$result_root/explicit/time.log")"
if [[ -z "$explicit_total" ]]; then
  echo "cannot parse explicit wall time" >&2
  exit 1
fi
explicit_peak="$(awk '/maximum resident set size/{print $1/1024/1024}' "$result_root/explicit/time.log" | tail -1)"
explicit_iter="$(python3 -c 'import pathlib,sys; rows=[x.split() for x in (pathlib.Path(sys.argv[1])/"Residuals.plt").read_text().splitlines() if x.strip() and not x.lstrip().startswith("VARIABLES")]; print(int(rows[-1][0]))' "$result_root/explicit")"
explicit_flow="$(python3 -c 'import pathlib,re,sys; files=list(pathlib.Path(sys.argv[1]).glob("Flow2d_[0-9]*.plt")); print(max(files,key=lambda p:int(re.search(r"_(\d+)\.plt$",p.name).group(1))))' "$result_root/explicit")"
printf 'Explicit,0,0,0,true,%s,%s,0,0,0,%s,%s,1\n' "$explicit_iter" "$explicit_iter" "$explicit_total" "$explicit_peak" >> "$csv"

for ncfl in 1 10 100 1000; do
  run_one ImplicitSteady "$ncfl" "implicit-$ncfl"
  run_dir="$result_root/implicit-$ncfl"
  metrics="$(python3 -c 'import math,pathlib,sys; rows=[x.split() for x in (pathlib.Path(sys.argv[1])/"ResidualImplicit.plt").read_text().splitlines() if x.strip() and not x.lstrip().startswith("VARIABLES")]; configured=float(sys.argv[2]); minimum=min(float(r[2]) for r in rows); retries=int(rows[-1][18]); (math.isclose(minimum,configured,rel_tol=1e-12,abs_tol=1e-14) and retries==0) or (_ for _ in ()).throw(SystemExit(f"pseudo-CFL fallback invalidates nCFL={configured:g}: minimum={minimum:g}, retries={retries}")); r=rows[-1]; print(",".join([f"{minimum:g}",str(retries),"true",r[0],r[0],r[14],r[15],r[16]]))' "$run_dir" "$ncfl")"
  implicit_flow="$(python3 -c 'import pathlib,re,sys; files=list(pathlib.Path(sys.argv[1]).glob("Flow2d_[0-9]*.plt")); print(max(files,key=lambda p:int(re.search(r"_(\d+)\.plt$",p.name).group(1))))' "$run_dir")"
  python3 "$repo_root/tests/helpers/compare_flow_fields.py" "$explicit_flow" "$implicit_flow" \
    --reference-surface "$result_root/explicit/Surf_CpCfCh.plt" \
    --candidate-surface "$run_dir/Surf_CpCfCh.plt"
  peak="$(awk '/maximum resident set size/{print $1/1024/1024}' "$run_dir/time.log" | tail -1)"
  total="$(awk '/ real /{print $1; exit}' "$run_dir/time.log")"
  if [[ -z "$total" ]]; then
    echo "cannot parse implicit wall time for nCFL=$ncfl" >&2
    exit 1
  fi
  speedup="$(python3 -c 'import sys; print(float(sys.argv[1])/float(sys.argv[2]))' "$explicit_total" "$total")"
  printf 'ImplicitSteady,%s,%s,%s,%s,%s\n' "$ncfl" "$metrics" "$total" "$peak" "$speedup" >> "$csv"
done
cat "$csv"
echo "results: $result_root"
