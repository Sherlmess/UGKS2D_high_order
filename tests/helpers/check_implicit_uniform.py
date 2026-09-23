#!/usr/bin/env python3
import math
import sys
from pathlib import Path


path = Path(sys.argv[1]) / "ResidualImplicit.plt"
if not path.exists():
    raise SystemExit("ResidualImplicit.plt was not created")
lines = path.read_text().splitlines()
header_lines = [line for line in lines if line.lstrip().startswith("VARIABLES")]
if len(header_lines) != 1:
    raise SystemExit(f"expected one implicit diagnostic header, got {len(header_lines)}")
columns = header_lines[0].split("=", 1)[1].split()
expected_columns = """iter dtPhysicalMin nCFL alpha rhoResidual rhoUResidual
rhoVResidual energyResidual hResidual bResidual momentError massDefect retries wallTime
prSeconds lusgsSeconds correctionExchangeSeconds evolutionSeconds retryTotal rhoRelative
rhoURelative rhoVRelative energyRelative rhoAbsolute rhoUAbsolute rhoVAbsolute energyAbsolute
hRelative bRelative hAbsolute bAbsolute microForwardSeconds microBackwardSeconds
microForwardBackwardExchangeSeconds""".split()
if columns != expected_columns:
    raise SystemExit(f"unexpected implicit diagnostic columns: {columns}")
rows = [line.split() for line in lines
        if line.strip() and not line.lstrip().startswith("VARIABLES")]
values = [float(value) for value in rows[-1]]
if len(values) != 34:
    raise SystemExit(f"expected 34 implicit diagnostic columns, got {len(values)}")
if not all(math.isfinite(value) for value in values):
    raise SystemExit("non-finite implicit residual record")
if max(values[4:11]) >= 1.0e-12:
    raise SystemExit(f"uniform residuals too large: {values}")
if int(values[18]) != 0:
    raise SystemExit(f"uniform state unexpectedly used pseudo-CFL retries: {values[18]}")
print("PASS implicit uniform zero-residual regression")
