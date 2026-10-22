#!/usr/bin/env bash
set -euo pipefail
code_dir="$(cd "$(dirname "$0")" && pwd)"
export S2_NEAR_WALL_VL=1 S2_NEAR_WALL_GRADED235=1
export R5_WALL_LIMIT="${R5_WALL_LIMIT:-2400}"
exec bash "$code_dir/run.sh" "$@"
