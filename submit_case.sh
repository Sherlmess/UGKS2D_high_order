#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
case_name="${1:?Use high_implicit, low_implicit or low_explicit}"
case "$case_name" in high_implicit|low_implicit|low_explicit) ;; *) exit 2;; esac
[[ -x "$root/bin/ugks_solver" ]] || { echo 'Compile in this NEW package first: bash build_hpc.sh' >&2; exit 2; }
command -v sbatch >/dev/null
mkdir -p "$root/runs"
run_dir=$(mktemp -d "$root/runs/${case_name}_$(date +%Y%m%dT%H%M%S)_XXXXXX")
cp "$root/cases/$case_name/Config.ini" "$root/cases/$case_name/Mesh2d.dat" "$root/cases/$case_name/velocity_half32_scale2.dat" "$run_dir/"
cp "$root/bin/ugks_solver" "$run_dir/ugks_solver"
chmod u+x "$run_dir/ugks_solver"
cp "$root/job.slurm" "$run_dir/job.slurm"
printf '%s\n' "$case_name" > "$run_dir/METHOD.txt"
cd "$run_dir"
sha256sum Config.ini Mesh2d.dat velocity_half32_scale2.dat ugks_solver job.slurm > INPUT_SHA256SUMS
printf 'Frozen run directory: %s\n' "$run_dir"
sbatch --job-name="C3_${case_name}" job.slurm | tee submission.txt
