#!/usr/bin/env python3
import json
import math
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
    data = values[-12 * ncell :]
    if len(data) != 12 * ncell:
        raise SystemExit(f"cannot parse {path}")
    return [data[k * ncell : (k + 1) * ncell] for k in range(12)]


def field_error(reference, candidate):
    # Use the four primitive disturbance blocks rho', u', v', and T'.
    diff2 = ref2 = 0.0
    linf = 0.0
    for k in range(2, 6):
        for ref, value in zip(reference[k], candidate[k]):
            delta = value - ref
            diff2 += delta * delta
            ref2 += ref * ref
            linf = max(linf, abs(delta))
    return math.sqrt(diff2 / max(ref2, 1.0e-60)), linf


root = Path(sys.argv[1])
reference = read_blocks(root / "reference" / "Flow2d_64.plt", 16, 4)
cases = [("dt0040", 2), ("dt0020", 4), ("dt0010", 8)]
errors = []
linfs = []
for name, step in cases:
    field = read_blocks(root / name / f"Flow2d_{step}.plt", 16, 4)
    error, linf = field_error(reference, field)
    errors.append(error)
    linfs.append(linf)

if not all(math.isfinite(value) and value > 0.0 for value in errors + linfs):
    raise SystemExit(f"nonfinite temporal errors: relative={errors}, linf={linfs}")
orders = [math.log(errors[i] / errors[i + 1], 2.0) for i in range(2)]
monotone = errors[0] > errors[1] > errors[2]
if not monotone:
    raise SystemExit(f"temporal errors are not monotone: {errors}")

# This is an evidence-producing refinement gate, not an unconditional order claim.
result = {
    "physical_dt": [0.04, 0.02, 0.01],
    "reference_dt": 0.00125,
    "final_time": 0.08,
    "relative_l2_rho_u_v_T": errors,
    "absolute_linf_rho_u_v_T": linfs,
    "observed_orders": orders,
    "formal_second_order_resolved": min(orders) >= 1.8,
}
(root / "temporal_metrics.json").write_text(json.dumps(result, indent=2) + "\n")
print("PASS implicit-unsteady temporal refinement")
print(json.dumps(result, indent=2))
