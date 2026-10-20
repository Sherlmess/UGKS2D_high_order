#!/usr/bin/env bash
#SBATCH --job-name=S2D_Implicit
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
# Site settings above are copied from the existing template; Intel/HPC is unverified here.
case_name="${CASE_NAME:-high_implicit}"
case "$case_name" in low_explicit|low_implicit|high_implicit) ;; *) exit 2 ;; esac
run_dir="run_${case_name}_${SLURM_JOB_ID}"
[[ ! -e "$run_dir" ]] || exit 2
mkdir "$run_dir"
cp "examples/$case_name/"{Config.ini,Mesh2d.dat,velocity_half32_scale2.dat} "$run_dir/"
bash run.sh "$run_dir"
