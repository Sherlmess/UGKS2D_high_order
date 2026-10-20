#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
method="${1:?Use high_implicit, low_implicit or low_explicit}"
case "$method" in high_implicit|low_implicit|low_explicit) ;; *) exit 2;; esac
[[ -x "$root/bin/ugks_solver" ]] || { echo 'Existing compiled bin/ugks_solver missing' >&2; exit 2; }
command -v sbatch >/dev/null
mkdir -p "$root/runs"
run=$(mktemp -d "$root/runs/smoke_${method}_XXXXXX")
cp "$root/cases/$method/"{Config.ini,Mesh2d.dat,velocity_half32_scale2.dat} "$run/"
# Two accepted steps maximum; unchanged physics, CFL, rank layout and scientific gates.
awk '/^[[:space:]]*(MAX_ITER|Implicit_MaxIter)[[:space:]]*=/ {sub(/=.*/,"= 2"); n++} {print} END{if(n!=2)exit 2}' "$run/Config.ini" > "$run/Config.ini.tmp"
mv "$run/Config.ini.tmp" "$run/Config.ini"
cp "$root/bin/ugks_solver" "$run/ugks_solver"
chmod u+x "$run/ugks_solver"
cp "$root/job.slurm" "$run/job.slurm"
printf '%s\n' "$method" > "$run/METHOD.txt"
cd "$run"
sha256sum Config.ini Mesh2d.dat velocity_half32_scale2.dat ugks_solver job.slurm > INPUT_SHA256SUMS
printf 'Two-step run directory: %s\n' "$run"
sbatch --time=00:15:00 --job-name="C3_smoke_${method}" job.slurm | tee submission.txt
