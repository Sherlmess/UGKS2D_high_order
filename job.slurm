#!/usr/bin/env bash
#SBATCH --partition=wzhcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.log
set -euo pipefail
module purge
module add compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -s unlimited
cd "${SLURM_SUBMIT_DIR:?}"
sha256sum -c INPUT_SHA256SUMS
get_int() { sed -nE "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*([0-9]+).*/\1/p" Config.ini; }
npx=$(get_int NP_X); npy=$(get_int NP_Y)
[[ "$npx" =~ ^[0-9]+$ && "$npy" =~ ^[0-9]+$ ]] || exit 2
np=$((npx*npy))
[[ "${SLURM_NTASKS:-0}" -eq "$np" ]] || { echo 'rank/config mismatch'; exit 2; }
# Prevent an inherited diagnostic/preflight environment from replacing formal work.
for key in $(compgen -e); do
 case "$key" in R4_*|R5_*|R6_*|W3B_*|NATIVE_*|S2_*) unset "$key";; esac
done
export S2_OWNER_RECON=1 S2_OWNER_COMPACT=1 S2_FIRST_LAYER_SCOPE=1 S2_VELOCITY_TILE=1 S2_WALL_TILE=1
# R6 soft stop applies to the closed implicit methods; original default was 230 s.
export R5_WALL_LIMIT=82800
printf 'method=%s ranks=%s start=%s\n' "$(cat METHOD.txt)" "$np" "$(date -u +%FT%TZ)" > job_status.txt
set +e
srun --mpi=pmi2 -n "$np" ./ugks_solver > stdout.log 2>&1
rc=$?
set -e
printf 'exit=%s finish=%s\n' "$rc" "$(date -u +%FT%TZ)" >> job_status.txt
shopt -s nullglob
refs=(reference_*.in)
if (( ${#refs[@]} )); then cp "${refs[${#refs[@]}-1]}" residual_reference.in; fi
exit "$rc"
