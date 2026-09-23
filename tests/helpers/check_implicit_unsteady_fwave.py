#!/usr/bin/env python3
import math
import re
import sys
from pathlib import Path


def read_blocks(path: Path, nx: int, ny: int):
    values = []
    for token in path.read_text().split():
        try:
            values.append(float(token.replace("D", "E").replace("d", "e")))
        except ValueError:
            pass
    ncell = nx * ny
    expected = 2 + 12 * ncell
    if len(values) < expected:
        raise SystemExit(f"not enough numeric data in {path}: {len(values)} < {expected}")
    # The first two numeric values are the I and J dimensions in the zone line.
    data = values[-12 * ncell :]
    return [data[k * ncell : (k + 1) * ncell] for k in range(12)]


def relative_l2(reference, candidate):
    numerator = math.sqrt(sum((a - b) ** 2 for a, b in zip(reference, candidate)))
    denominator = math.sqrt(sum(a * a for a in reference))
    return numerator / max(denominator, 1.0e-30)


def absolute_linf(reference, candidate):
    return max(abs(a - b) for a, b in zip(reference, candidate))


root = Path(sys.argv[1])
explicit_dir = root / "explicit"
implicit_dir = root / "implicit"
explicit = read_blocks(explicit_dir / "Flow2d_1.plt", 16, 4)
implicit = read_blocks(implicit_dir / "Flow2d_1.plt", 16, 4)

# Blocks 2--5 are rho', u', v', and T' after x and y.
errors = [relative_l2(explicit[k], implicit[k]) for k in range(2, 6)]
absolute_errors = [absolute_linf(explicit[k], implicit[k]) for k in range(2, 6)]
signal = max(abs(value) for block in explicit[2:6] for value in block)
if not math.isfinite(signal) or signal < 1.0e-7:
    raise SystemExit(f"F-wave response is missing or nonfinite: max signal={signal:.6e}")
if not all(math.isfinite(value) and value <= 3.0e-3 for value in errors):
    raise SystemExit(f"implicit explicit-limit relative mismatch for rho/u/v/T: {errors}")
if not all(math.isfinite(value) and value <= 5.0e-8 for value in absolute_errors):
    raise SystemExit(f"implicit explicit-limit absolute mismatch for rho/u/v/T: {absolute_errors}")

log = (implicit_dir / "solver.log").read_text()
if "Unsteady implicit step rejected" in log or "ERROR" in log:
    raise SystemExit("implicit run rejected a physical step")
rows = re.findall(r"Unsteady step=(\d+) inner=(\d+)", log)
if len(rows) != 1 or max(int(inner) for _, inner in rows) > 1:
    raise SystemExit(f"unexpected implicit step diagnostics: {rows}")
print("PASS implicit-unsteady F-wave explicit limit")
print("relative L2 rho/u/v/T:", " ".join(f"{value:.6e}" for value in errors))
print("absolute Linf rho/u/v/T:", " ".join(f"{value:.6e}" for value in absolute_errors))
print("inner iterations:", rows[0][1])
