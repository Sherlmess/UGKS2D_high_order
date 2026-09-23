#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
grid="${TEMPORAL_GRID:-96}"

case "$grid" in
    96)
        nx=96
        ny=6
        reference_dt=0.0003125
        reference_steps=128
        refinement=0
        ;;
    160)
        nx=160
        ny=8
        reference_dt=0.00015625
        reference_steps=256
        refinement=1
        ;;
    *)
        echo "TEMPORAL_GRID must be 96 or 160" >&2
        exit 2
        ;;
esac

if [[ -n "${TEMPORAL_RUN_ROOT:-}" ]]; then
    run_root="$TEMPORAL_RUN_ROOT"
    mkdir -p "$run_root"
else
    run_root="$(mktemp -d /private/tmp/implicit-unsteady-temporal-smooth.XXXXXX)"
fi
if find "$run_root" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
    echo "refusing to reuse non-empty temporal output directory: $run_root" >&2
    exit 2
fi

mkdir -p "$run_root/base" "$run_root/reference" \
    "$run_root/dt0010" "$run_root/dt0005" "$run_root/dt00025" "$run_root/dt000125"
printf 'case\truntime_seconds\n' > "$run_root/timings.tsv"
experiment_start="$(python3 -c 'import time; print(time.monotonic())')"

run_case() {
    local case_dir="$1"
    if ! python3 - "$case_dir" "$run_root" "$experiment_start" <<'PY'
import pathlib
import subprocess
import sys
import time

case_dir = pathlib.Path(sys.argv[1])
run_root = pathlib.Path(sys.argv[2])
experiment_start = float(sys.argv[3])
deadline = 30.0 * 60.0
elapsed_before = time.monotonic() - experiment_start
if elapsed_before >= deadline:
    print("temporal experiment exceeded 30 minutes before " + case_dir.name, file=sys.stderr)
    raise SystemExit(124)

started = time.monotonic()
try:
    with (case_dir / "solver.log").open("w") as log:
        result = subprocess.run(
            ["mpirun", "-np", "1", "./ugks_solver"],
            cwd=case_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=max(1.0, deadline - elapsed_before),
            check=False,
        )
except subprocess.TimeoutExpired:
    print("temporal experiment exceeded 30 minutes in " + case_dir.name, file=sys.stderr)
    raise SystemExit(124)

runtime = time.monotonic() - started
(case_dir / "runtime_seconds.txt").write_text(f"{runtime:.9f}\n")
with (run_root / "timings.tsv").open("a") as timing:
    timing.write(f"{case_dir.name}\t{runtime:.9f}\n")
if result.returncode != 0:
    print(f"{case_dir.name} exited with {result.returncode}", file=sys.stderr)
    raise SystemExit(result.returncode)
if time.monotonic() - experiment_start > deadline:
    print("temporal experiment exceeded 30 minutes after " + case_dir.name, file=sys.stderr)
    raise SystemExit(124)
PY
    then
        cat "$case_dir/solver.log" >&2
        echo "temporal run root: $run_root" >&2
        exit 1
    fi
}

cp "$repo_root/bin/ugks_solver" \
    "$repo_root/tests/cases/implicit_unsteady_uniform/Config.ini" "$run_root/base/"
python3 "$repo_root/tests/helpers/make_cartesian_mesh.py" "$run_root/base/Mesh2d.dat" \
    --nx "$nx" --ny "$ny" --surf "$run_root/base/surf.dat"
python3 "$repo_root/tests/helpers/patch_config.py" "$run_root/base/Config.ini" \
    "TimeIntegration_method='Explicit'" Iflag_InitFlow=0 Iflag_WriteSurf=0 Iflag_WriteRestart=1 \
    Nx_all="$nx" Ny_all="$ny" Nu=16 Nv=16 NP_X=1 NP_Y=1 Quadrature_method=2 \
    Ma=0.1d0 Re=100.0d0 Tinf=1.0d0 Tw=1.0d0 Pr=1.0d0 \
    BC_Type="4 4 6 6" FluxWall_type=0 ACC=1.0d0 TimeStep_method=0 \
    dt_global="${reference_dt}d0" CFL=0.025d0 MAX_ITER=1 Nstep_show=1 Nstep_write=1
run_case "$run_root/base"

cp "$run_root/base/flow2d.dat" "$run_root/uniform_flow2d.dat"
cp "$run_root/base/vdf2d.dat" "$run_root/uniform_vdf2d.dat"
cp "$run_root/base/Mesh2d.dat" "$run_root/initial_Mesh2d.dat"
python3 "$repo_root/tests/helpers/make_smooth_bump_restart.py" \
    --mesh "$run_root/base/Mesh2d.dat" \
    --flow-in "$run_root/base/flow2d.dat" \
    --vdf-in "$run_root/base/vdf2d.dat" \
    --flow-out "$run_root/smooth_flow2d.dat" \
    --vdf-out "$run_root/smooth_vdf2d.dat" \
    --metadata "$run_root/smooth_bump_metadata.json"
cp "$run_root/smooth_flow2d.dat" "$run_root/initial_flow2d.dat"
cp "$run_root/smooth_vdf2d.dat" "$run_root/initial_vdf2d.dat"

python3 - "$run_root/temporal_matrix.json" "$nx" "$ny" "$reference_dt" "$reference_steps" "$refinement" <<'PY'
import json
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
matrix = {
    "experiment": "smooth_compatible_bump",
    "nx": int(sys.argv[2]),
    "ny": int(sys.argv[3]),
    "candidate_dt": [0.01, 0.005, 0.0025, 0.00125],
    "candidate_steps": [4, 8, 16, 32],
    "reference_dt": float(sys.argv[4]),
    "reference_steps": int(sys.argv[5]),
    "reference_method": "ImplicitUnsteady",
    "final_time": 0.04,
    "theta": 0.5,
    "inner_tolerance": 1.0e-10,
    "max_inner": 100,
    "cfl": 0.025,
    "iflag_dis": 0,
    "bc_type": [4, 4, 6, 6],
    "quadrature_method": 2,
    "nu": 16,
    "nv": 16,
    "perturbation_interval": [0.25, 0.75],
    "perturbation_amplitude": 1.0e-5,
    "refinement": bool(int(sys.argv[6])),
}
path.write_text(json.dumps(matrix, indent=2) + "\n")
PY

make_case() {
    local name="$1"
    local dt="$2"
    local steps="$3"
    local case_dir="$run_root/$name"

    cp "$repo_root/bin/ugks_solver" "$repo_root/tests/cases/implicit_unsteady_uniform/Config.ini" \
        "$run_root/base/Mesh2d.dat" "$run_root/base/surf.dat" \
        "$run_root/smooth_flow2d.dat" "$run_root/smooth_vdf2d.dat" "$case_dir/"
    mv "$case_dir/smooth_flow2d.dat" "$case_dir/flow2d.dat"
    mv "$case_dir/smooth_vdf2d.dat" "$case_dir/vdf2d.dat"
    cp "$case_dir/flow2d.dat" "$case_dir/input_flow2d.dat"
    cp "$case_dir/vdf2d.dat" "$case_dir/input_vdf2d.dat"
    python3 "$repo_root/tests/helpers/patch_config.py" "$case_dir/Config.ini" \
        "TimeIntegration_method='ImplicitUnsteady'" Iflag_InitFlow=2 Iflag_WriteSurf=0 Iflag_WriteRestart=0 \
        Nx_all="$nx" Ny_all="$ny" Nu=16 Nv=16 NP_X=1 NP_Y=1 Quadrature_method=2 \
        Ma=0.1d0 Re=100.0d0 Tinf=1.0d0 Tw=1.0d0 Pr=1.0d0 \
        BC_Type="4 4 6 6" FluxWall_type=0 ACC=1.0d0 TimeStep_method=0 \
        dt_global="${dt}d0" CFL=0.025d0 MAX_ITER="$steps" Nstep_show=1 Nstep_write="$steps" \
        Unsteady_Theta=0.5d0 Unsteady_InnerTol=1.0d-10 Unsteady_MaxInner=100
}

make_case reference "$reference_dt" "$reference_steps"
make_case dt0010 0.01 4
make_case dt0005 0.005 8
make_case dt00025 0.0025 16
make_case dt000125 0.00125 32

for name in reference dt0010 dt0005 dt00025 dt000125; do
    run_case "$run_root/$name"
done

checker_rc=0
python3 "$repo_root/tests/helpers/check_implicit_unsteady_temporal_order.py" "$run_root" || checker_rc=$?
echo "temporal run root: $run_root"
exit "$checker_rc"
