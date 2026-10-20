#!/usr/bin/env bash
set -euo pipefail
code_dir="$(cd "$(dirname "$0")" && pwd)"
run_dir="${1:?Usage: bash run.sh PREPARED_RUN_DIRECTORY}"
[[ -f "$run_dir/Config.ini" && -f "$run_dir/Mesh2d.dat" ]] || { echo 'Prepare Config.ini and mesh first' >&2; exit 2; }
[[ -x "$code_dir/bin/ugks_solver" ]] || { echo 'Run build.sh first' >&2; exit 2; }
get_int() { sed -nE "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*([0-9]+).*/\1/p" "$run_dir/Config.ini"; }
npx=$(get_int NP_X); npy=$(get_int NP_Y)
[[ "$npx" =~ ^[0-9]+$ && "$npy" =~ ^[0-9]+$ ]] || exit 2
np=$((npx*npy))
[[ -z "${SLURM_NTASKS:-}" || "$SLURM_NTASKS" -eq "$np" ]] || { echo 'SLURM_NTASKS must equal NP_X*NP_Y' >&2; exit 2; }
cp "$code_dir/bin/ugks_solver" "$run_dir/ugks_solver"
cd "$run_dir"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# O13 switches select the already verified implementation; numerical formulas are unchanged.
export S2_OWNER_RECON=1 S2_OWNER_COMPACT=1 S2_FIRST_LAYER_SCOPE=1 S2_VELOCITY_TILE=1 S2_WALL_TILE=1
if [[ "${MPI_OVERSUBSCRIBE:-0}" == 1 ]]; then
 mpirun --oversubscribe -np "$np" ./ugks_solver > stdout.log 2>&1
else
 mpirun -np "$np" ./ugks_solver > stdout.log 2>&1
fi
# Carry the established micro residual reference epoch into the next complete restart.
shopt -s nullglob
refs=(reference_*.in)
if [[ -n "${refs[0]:-}" ]]; then cp "${refs[${#refs[@]}-1]}" residual_reference.in; fi
