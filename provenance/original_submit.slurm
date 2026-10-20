#!/usr/bin/env bash
#SBATCH --job-name=VL_Couette_orig
#SBATCH --partition=wzhcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=couette-original-%j.log
set -euo pipefail
module purge
module add compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -s unlimited
cd "${SLURM_SUBMIT_DIR:?Submit from package directory}"
# Each job copies its executable and input into a new directory before running.
bash run.sh "run_original_direction_${SLURM_JOB_ID}"
