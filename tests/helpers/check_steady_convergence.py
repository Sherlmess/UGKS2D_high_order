#!/usr/bin/env python3
import argparse
import math
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=["Explicit", "ImplicitSteady"])
parser.add_argument("run_dir", type=Path)
parser.add_argument("--tol", type=float, default=1.0e-8)
args = parser.parse_args()
path = args.run_dir / ("Residuals.plt" if args.mode == "Explicit" else "ResidualImplicit.plt")
status_path = args.run_dir / "ConvergenceStatus.txt"
if not status_path.exists():
    raise SystemExit(f"missing convergence status: {status_path}")
status_rows = [line.split() for line in status_path.read_text().splitlines() if line.strip()]
if len(status_rows) < 2 or status_rows[-1][1].lower() != "t":
    raise SystemExit(f"{args.mode} terminated without convergence: {status_rows[-1] if status_rows else []}")
if not path.exists():
    raise SystemExit(f"missing convergence file: {path}")
rows = [line.split() for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("VARIABLES")]
if not rows:
    raise SystemExit(f"empty convergence file: {path}")
values = [float(value) for value in rows[-1]]
residuals = values[1:5] if args.mode == "Explicit" else values[4:11]
if not all(math.isfinite(value) for value in residuals) or max(residuals) >= args.tol:
    raise SystemExit(f"{args.mode} did not converge: {residuals}")
print(f"PASS {args.mode} steady convergence")
