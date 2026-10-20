#!/usr/bin/env bash
set -euo pipefail
code_dir="$(cd "$(dirname "$0")" && pwd)"
old_dir="${1:?Usage: bash resume.sh OLD_COMPLETE_RUN NEW_DIRECTORY [ADDITIONAL_STEPS]}"
new_dir="${2:?New directory required}"
steps="${3:-3}"
[[ "$steps" =~ ^[1-9][0-9]*$ ]] || exit 2
[[ ! -e "$new_dir" ]] || { echo 'Refusing to overwrite an existing run' >&2; exit 2; }
for file in Config.ini Mesh2d.dat velocity_half32_scale2.dat flow2d.dat vdf2d.dat residual_reference.in; do
 [[ -f "$old_dir/$file" ]] || { echo "Missing complete restart file: $file" >&2; exit 2; }
done
mkdir "$new_dir"
cp "$old_dir"/{Config.ini,Mesh2d.dat,velocity_half32_scale2.dat,flow2d.dat,vdf2d.dat,residual_reference.in} "$new_dir/"
shopt -s nullglob
companions=("$old_dir"/closed_mass_rank*.bin)
[[ -n "${companions[0]:-}" ]] || { echo 'Missing R6 mass companions' >&2; exit 2; }
cp "${companions[@]}" "$new_dir/"
sed -E -e 's/^(Iflag_InitFlow[[:space:]]*=[[:space:]]*)[0-9]+/\12/' \
 -e "s/^(Implicit_MaxIter[[:space:]]*=[[:space:]]*)[0-9]+/\1$steps/" \
 -e "s/^(MAX_ITER[[:space:]]*=[[:space:]]*)[0-9]+/\1$steps/" \
 "$new_dir/Config.ini" > "$new_dir/Config.ini.next"
mv "$new_dir/Config.ini.next" "$new_dir/Config.ini"
bash "$code_dir/run.sh" "$new_dir"
