#!/usr/bin/env python3
import math
import re
import sys
from pathlib import Path

run_dir = Path(sys.argv[1])
log = (run_dir / 'solver.log').read_text()
rows = re.findall(
    r'Unsteady step=(\d+) inner=(\d+) Rrho/RrhoU/RrhoV/RE/Rh/Rb=\s*'
    r'([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+'
    r'([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)', log)
if len(rows) != 2:
    raise SystemExit(f'expected two accepted physical-step diagnostics, got {len(rows)}')
for row in rows:
    values = [float(value) for value in row[2:]]
    if not all(math.isfinite(value) and value <= 1.0e-8 for value in values):
        raise SystemExit(f'nonzero/nonfinite uniform residuals: {values}')
if int(rows[0][1]) > 1 or int(rows[1][1]) != 0:
    raise SystemExit(f'uniform state should require at most one startup equilibration correction: {rows}')
status = (run_dir / 'ConvergenceStatus.txt').read_text()
if 'ImplicitUnsteady F 2' not in status:
    raise SystemExit(f'unexpected physical-step termination status: {status}')
flow_files = sorted(run_dir.glob('Flow2d_*.plt'))
if not flow_files or not all(math.isfinite(float(token.replace('D','E')))
    for token in flow_files[-1].read_text().split()
    if re.fullmatch(r'[+\-0-9.EeDd]+', token)):
    raise SystemExit('missing or nonfinite final uniform field')
print('PASS implicit unsteady GH16 uniform regression')
