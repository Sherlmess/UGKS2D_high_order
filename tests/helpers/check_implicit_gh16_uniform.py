#!/usr/bin/env python3
import math
import sys
from pathlib import Path


run_dir = Path(sys.argv[1])
status_path = run_dir / "ConvergenceStatus.txt"
residual_path = run_dir / "ResidualImplicit.plt"
log_path = run_dir / "solver.log"
if not status_path.exists() or not residual_path.exists():
    raise SystemExit("missing implicit convergence output")
if not log_path.exists():
    raise SystemExit("missing captured implicit solver output")
log_text = log_path.read_text()
required_tokens = (
    "Implicit iter=", " nCFL=", " alpha=", " Rrho=", " RrhoU=",
    " RrhoV=", " RE=", " Rh=", " Rb=", " moment=", " mass=",
    " wall=", " retries=", " retryTotal=",
)
if not all(token in log_text for token in required_tokens):
    raise SystemExit("missing complete implicit console residual record")

status_rows = [line.split() for line in status_path.read_text().splitlines() if line.strip()]
if len(status_rows) < 2:
    raise SystemExit("missing implicit convergence status row")
status = status_rows[-1]
if len(status) < 4 or status[0] != "ImplicitSteady" or status[1].lower() != "t":
    raise SystemExit(f"GH16 uniform flow did not converge: {status}")
iteration = int(status[2])
if iteration < 1 or iteration > 5:
    raise SystemExit(f"GH16 uniform flow exceeded five iterations: {iteration}")

rows = [line.split() for line in residual_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("VARIABLES")]
values = [float(value) for value in rows[-1]]
if len(values) != 31 or not all(math.isfinite(value) for value in values):
    raise SystemExit("GH16 uniform flow has an invalid implicit residual record")
if max(values[4:11]) >= 1.0e-12:
    raise SystemExit(f"GH16 uniform residuals are too large: {values}")
if int(values[18]) != 0:
    raise SystemExit(f"GH16 uniform flow used pseudo-CFL retries: {values[18]}")

flow_path = run_dir / f"Flow2d_{iteration}.plt"
if not flow_path.exists():
    raise SystemExit(f"missing terminal flow output: {flow_path.name}")
numeric_values = []
for line in flow_path.read_text().splitlines():
    stripped = line.strip()
    if not stripped or stripped.startswith("VARIABLES") or stripped.startswith("ZONE"):
        continue
    try:
        numeric_values.extend(float(token) for token in stripped.split())
    except ValueError as exc:
        raise SystemExit(f"invalid terminal flow record: {stripped}") from exc
if not numeric_values or not all(math.isfinite(value) for value in numeric_values):
    raise SystemExit("GH16 uniform flow contains non-finite terminal field values")

print("PASS implicit GH16 uniform zero-residual regression")
