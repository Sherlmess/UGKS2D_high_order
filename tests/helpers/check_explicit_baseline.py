#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("run_dir", type=Path)
parser.add_argument("baseline", type=Path)
parser.add_argument("--record", action="store_true")
args = parser.parse_args()
flow = args.run_dir / "Flow2d_3.plt"
residual = args.run_dir / "Residuals.plt"
if not flow.exists() or not residual.exists():
    raise SystemExit("missing Flow2d_3.plt or Residuals.plt")
sha = hashlib.sha256(flow.read_bytes()).hexdigest()
rows = [
    line.split()
    for line in residual.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("VARIABLES")
]
values = [float(value) for value in rows[-1][1:5]]
record = {"flow_sha256": sha, "iteration3_residuals": values}
if args.record:
    args.baseline.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
else:
    expected = json.loads(args.baseline.read_text())
    if record != expected:
        raise SystemExit(
            f"explicit baseline mismatch\nexpected={expected}\nactual={record}"
        )
print("PASS explicit structured three-step baseline")
