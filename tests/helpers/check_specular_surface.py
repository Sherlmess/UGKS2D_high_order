#!/usr/bin/env python3
import argparse
import math
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("surface", type=Path)
parser.add_argument("--tol", type=float, default=1.0e-8)
args = parser.parse_args()
rows = [[float(token) for token in line.replace("D", "E").split()]
        for line in args.surface.read_text().splitlines()
        if line.strip() and not line.lstrip().lower().startswith("variables")]
if not rows or not all(len(row) == 5 and all(math.isfinite(value) for value in row) for row in rows):
    raise SystemExit(f"invalid surface data: {args.surface}")
max_cf = max(abs(row[3]) for row in rows)
max_ch = max(abs(row[4]) for row in rows)
print(f"specular surface max|Cf|={max_cf:.8e} max|Ch|={max_ch:.8e}")
if max_cf >= args.tol or max_ch >= args.tol:
    raise SystemExit("stationary specular wall has nonzero shear or heat transfer")
print("PASS stationary specular surface invariants")
