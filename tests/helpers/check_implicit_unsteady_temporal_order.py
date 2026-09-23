#!/usr/bin/env python3
import json
import math
import re
import struct
import sys
from array import array
from pathlib import Path


GH16_COORDS = [
    -0.3686007162724397e1, -0.2863133883708075e1,
    -0.2183921153095858e1, -0.1588855862270055e1,
    -0.1064246312116224e1, -0.6163028841823999,
    -0.2673983721677653, -0.5297864393185113e-1,
    0.5297864393185113e-1, 0.2673983721677653,
    0.6163028841823999, 0.1064246312116224e1,
    0.1588855862270055e1, 0.2183921153095858e1,
    0.2863133883708075e1, 0.3686007162724397e1,
]
GH16_WEIGHTS = [
    0.1192596926595344e-5, 0.2020636491324107e-3,
    0.5367935756025333e-2, 0.4481410991746290e-1,
    0.1574482826187903, 0.2759533979884218,
    0.2683307544726388, 0.1341091884533595,
    0.1341091884533595, 0.2683307544726388,
    0.2759533979884218, 0.1574482826187903,
    0.4481410991746290e-1, 0.5367935756025333e-2,
    0.2020636491324107e-3, 0.1192596926595344e-5,
]


def fail(message, code=2):
    print(message, file=sys.stderr)
    raise SystemExit(code)


def as_float(value):
    return float(value.replace("D", "E").replace("d", "e"))


def config_value(path, key):
    pattern = re.compile(rf"(?mi)^\s*{re.escape(key)}\s*=\s*([^!\r\n]+)")
    match = pattern.search(path.read_text())
    if not match:
        fail(f"missing {key} in {path}")
    return match.group(1).strip().split()[0].rstrip(",/")


def config_real(path, key):
    return as_float(config_value(path, key))


def config_int(path, key):
    return int(config_value(path, key))


def read_records(path):
    data = path.read_bytes()
    records = []
    offset = 0
    while offset < len(data):
        if offset + 4 > len(data):
            fail(f"truncated record marker in {path}")
        length = struct.unpack_from("=i", data, offset)[0]
        if length < 0:
            fail(f"negative record length in {path}")
        start = offset + 4
        end = start + length
        if end + 4 > len(data):
            fail(f"truncated record payload in {path}")
        if struct.unpack_from("=i", data, end)[0] != length:
            fail(f"record marker mismatch in {path}")
        records.append(data[start:end])
        offset = end + 4
    return records


def read_ints(payload, count, label):
    if len(payload) != 4 * count:
        fail(f"{label} has {len(payload)} bytes, expected {4 * count}")
    return struct.unpack("=" + "i" * count, payload)


def read_doubles(payload, count, label):
    if len(payload) != 8 * count:
        fail(f"{label} has {len(payload)} bytes, expected {8 * count}")
    values = array("d")
    values.frombytes(payload)
    if len(values) != count:
        fail(f"{label} has {len(values)} doubles, expected {count}")
    return list(values)


def read_mesh(path, nx, ny):
    records = read_records(path)
    if len(records) != 2:
        fail(f"{path} must contain two records")
    nxp, nyp = read_ints(records[0], 2, "mesh dimensions")
    if (nxp, nyp) != (nx + 1, ny + 1):
        fail(f"mesh dimensions are {(nxp - 1, nyp - 1)}, expected {(nx, ny)}")
    node_count = nxp * nyp
    coords = read_doubles(records[1], 2 * node_count, "mesh coordinates")
    x_nodes = coords[:node_count]
    y_nodes = coords[node_count:]
    cell_x = []
    cell_y = []
    areas = []
    length1 = []
    length2 = []
    for j in range(ny):
        for i in range(nx):
            lower_left = j * nxp + i
            lower_right = lower_left + 1
            upper_left = (j + 1) * nxp + i
            upper_right = upper_left + 1
            cell_x.append(0.25 * (x_nodes[lower_left] + x_nodes[lower_right] +
                                  x_nodes[upper_left] + x_nodes[upper_right]))
            cell_y.append(0.25 * (y_nodes[lower_left] + y_nodes[lower_right] +
                                  y_nodes[upper_left] + y_nodes[upper_right]))
            bottom = math.hypot(x_nodes[lower_right] - x_nodes[lower_left],
                                y_nodes[lower_right] - y_nodes[lower_left])
            top = math.hypot(x_nodes[upper_right] - x_nodes[upper_left],
                             y_nodes[upper_right] - y_nodes[upper_left])
            left = math.hypot(x_nodes[upper_left] - x_nodes[lower_left],
                              y_nodes[upper_left] - y_nodes[lower_left])
            right = math.hypot(x_nodes[upper_right] - x_nodes[lower_right],
                               y_nodes[upper_right] - y_nodes[lower_right])
            length1.append(0.5 * (bottom + top))
            length2.append(0.5 * (left + right))
            areas.append(abs((x_nodes[lower_left] - x_nodes[upper_right]) *
                             (y_nodes[lower_right] - y_nodes[upper_left]) -
                             (x_nodes[lower_right] - x_nodes[upper_left]) *
                             (y_nodes[lower_left] - y_nodes[upper_right])) * 0.5)
    return cell_x, cell_y, areas, length1, length2


def read_restart(flow_path, vdf_path, nx, ny):
    cells = nx * ny
    flow_records = read_records(flow_path)
    if len(flow_records) != 3:
        fail(f"{flow_path} must contain three records")
    if read_ints(flow_records[0], 1, "flow header") != (1,):
        fail(f"unexpected flow header in {flow_path}")
    flow_nx, flow_ny, nvar = read_ints(flow_records[1], 3, "flow dimensions")
    if (flow_nx, flow_ny, nvar) != (nx, ny, 4):
        fail(f"unexpected flow dimensions in {flow_path}: {(flow_nx, flow_ny, nvar)}")
    flow = read_doubles(flow_records[2], cells * 4, "flow data")

    vdf_records = read_records(vdf_path)
    if len(vdf_records) != 3:
        fail(f"{vdf_path} must contain three records")
    if read_ints(vdf_records[0], 1, "vdf header") != (1,):
        fail(f"unexpected vdf header in {vdf_path}")
    vdf_nx, vdf_ny, vdf_nvar = read_ints(vdf_records[1], 3, "vdf dimensions")
    if (vdf_nx, vdf_ny) != (nx, ny) or vdf_nvar % 2 != 0:
        fail(f"unexpected vdf dimensions in {vdf_path}: {(vdf_nx, vdf_ny, vdf_nvar)}")
    nvel = vdf_nvar // 2
    if nvel != 256:
        fail(f"expected GH16 VDF with 256 velocities, found {nvel}")
    vdf = read_doubles(vdf_records[2], cells * vdf_nvar, "vdf data")
    return flow, vdf, nvel


def gh16_velocity_space():
    u = []
    v = []
    weight = []
    for j, vcoord in enumerate(GH16_COORDS):
        for i, ucoord in enumerate(GH16_COORDS):
            u.append(ucoord)
            v.append(vcoord)
            weight.append(
                GH16_WEIGHTS[i] * math.exp(ucoord * ucoord) *
                GH16_WEIGHTS[j] * math.exp(vcoord * vcoord)
            )
    return u, v, weight


def moment_values(h, b, u, v, weight):
    return [
        sum(weight[k] * h[k] for k in range(len(h))),
        sum(weight[k] * u[k] * h[k] for k in range(len(h))),
        sum(weight[k] * v[k] * h[k] for k in range(len(h))),
        0.5 * sum(weight[k] * ((u[k] * u[k] + v[k] * v[k]) * h[k] + b[k])
                  for k in range(len(h))),
    ]


def verify_restart_compatibility(flow, vdf, nx, ny, nvel):
    cells = nx * ny
    u, v, weight = gh16_velocity_space()
    max_absolute = 0.0
    max_relative = 0.0
    for cell in range(cells):
        target = [flow[component * cells + cell] for component in range(4)]
        h = [vdf[velocity * cells + cell] for velocity in range(nvel)]
        b = [vdf[(nvel + velocity) * cells + cell] for velocity in range(nvel)]
        if not all(math.isfinite(value) for value in target + h + b):
            fail("nonfinite smooth bump restart value")
        computed = moment_values(h, b, u, v, weight)
        scale = max(1.0, max(abs(value) for value in target))
        for actual, expected in zip(computed, target):
            error = abs(actual - expected)
            max_absolute = max(max_absolute, error)
            max_relative = max(max_relative, error / scale)
    tolerance = 1.0e-10
    if max_relative > tolerance:
        fail(f"smooth bump restart moment mismatch: {max_relative:.6e} > {tolerance:.6e}")
    return {
        "max_absolute_error": max_absolute,
        "max_relative_error": max_relative,
        "tolerance": tolerance,
        "passed": True,
    }


def read_plt(path, nx, ny):
    values = []
    for token in path.read_text().split():
        try:
            values.append(as_float(token))
        except ValueError:
            pass
    ncell = nx * ny
    expected = 12 * ncell
    if len(values) < expected:
        fail(f"not enough numeric data in {path}: {len(values)} < {expected}")
    data = values[-expected:]
    fields = [data[index * ncell:(index + 1) * ncell] for index in range(12)]
    if not all(math.isfinite(value) for field in fields for value in field):
        fail(f"nonfinite field in {path}")
    return fields


def primitive_from_flow(flow, nx, ny):
    cells = nx * ny
    gamma = 5.0 / 3.0
    result = [[] for _ in range(4)]
    for cell in range(cells):
        rho = flow[cell]
        rho_u = flow[cells + cell]
        rho_v = flow[2 * cells + cell]
        energy = flow[3 * cells + cell]
        if not all(math.isfinite(value) for value in (rho, rho_u, rho_v, energy)) or rho <= 0.0:
            fail("invalid initial smooth bump conservative state")
        u = rho_u / rho
        v = rho_v / rho
        denominator = energy - 0.5 * (rho_u * rho_u + rho_v * rho_v) / rho
        if denominator <= 0.0 or not math.isfinite(denominator):
            fail("invalid initial smooth bump thermal energy")
        inverse_temperature = 0.5 * rho / ((gamma - 1.0) * denominator)
        temperature = 1.0 / inverse_temperature
        result[0].append(rho)
        result[1].append(u)
        result[2].append(v)
        result[3].append(temperature)
    return result


def field_metric(reference, candidate):
    diff2 = 0.0
    reference2 = 0.0
    absolute = 0.0
    for expected, actual in zip(reference, candidate):
        if not math.isfinite(expected) or not math.isfinite(actual):
            fail("nonfinite temporal perturbation field")
        difference = actual - expected
        diff2 += difference * difference
        reference2 += expected * expected
        absolute = max(absolute, abs(difference))
    reference_norm = math.sqrt(reference2)
    return {
        "relative_l2": math.sqrt(diff2) / max(reference_norm, 1.0e-300),
        "absolute_linf": absolute,
        "reference_norm": reference_norm,
    }


def combined_metric(reference, candidate):
    diff2 = 0.0
    reference2 = 0.0
    absolute = 0.0
    for reference_field, candidate_field in zip(reference, candidate):
        for expected, actual in zip(reference_field, candidate_field):
            if not math.isfinite(expected) or not math.isfinite(actual):
                fail("nonfinite temporal perturbation field")
            difference = actual - expected
            diff2 += difference * difference
            reference2 += expected * expected
            absolute = max(absolute, abs(difference))
    reference_norm = math.sqrt(reference2)
    return {
        "relative_l2": math.sqrt(diff2) / max(reference_norm, 1.0e-300),
        "absolute_linf": absolute,
        "reference_norm": reference_norm,
    }


def observed_orders(errors):
    result = []
    smallest = sys.float_info.min
    for left, right in zip(errors, errors[1:]):
        if not math.isfinite(left) or not math.isfinite(right) or left < 0.0 or right < 0.0:
            fail(f"invalid error sequence: {errors}")
        result.append(math.log(max(left, smallest) / max(right, smallest), 2.0))
    return result


def check_dis_controls(path):
    lines = path.read_text().splitlines()
    try:
        marker = next(index for index, line in enumerate(lines)
                      if line.strip() == "$DIS_para")
    except StopIteration:
        fail(f"missing $DIS_para in {path}")
    if marker + 3 >= len(lines):
        fail(f"truncated DIS controls in {path}")
    values = [int(lines[marker + offset].strip()) for offset in range(1, 4)]
    if values != [0, 0, 0]:
        fail(f"unexpected DIS controls in {path}: {values}")
    if any(line.strip().startswith("F ") for line in lines):
        fail(f"sudden F-wave remains in {path}")


def check_config(path, matrix, dt, steps):
    expected_ints = {
        "Nx_all": matrix["nx"],
        "Ny_all": matrix["ny"],
        "Nu": 16,
        "Nv": 16,
        "NP_X": 1,
        "NP_Y": 1,
        "Quadrature_method": 2,
        "TimeStep_method": 0,
        "Iflag_InitFlow": 2,
        "Iflag_WriteSurf": 0,
        "Iflag_WriteRestart": 0,
        "MAX_ITER": steps,
        "Nstep_show": 1,
        "Nstep_write": steps,
        "FluxWall_type": 0,
        "Unsteady_MaxInner": 100,
    }
    for key, expected in expected_ints.items():
        if config_int(path, key) != expected:
            fail(f"{path} has unexpected {key}")
    if config_value(path, "TimeIntegration_method").strip("'\"") != "ImplicitUnsteady":
        fail(f"{path} is not ImplicitUnsteady")
    if not math.isclose(config_real(path, "dt_global"), dt, rel_tol=0.0, abs_tol=1.0e-15):
        fail(f"{path} has unexpected dt_global")
    expected_reals = {
        "Ma": 0.1,
        "Re": 100.0,
        "Tinf": 1.0,
        "Tw": 1.0,
        "Pr": 1.0,
        "CFL": matrix["cfl"],
        "Unsteady_Theta": matrix["theta"],
        "Unsteady_InnerTol": matrix["inner_tolerance"],
    }
    for key, expected in expected_reals.items():
        if not math.isclose(config_real(path, key), expected, rel_tol=0.0, abs_tol=1.0e-14):
            fail(f"{path} has unexpected {key}")
    text = path.read_text()
    match = re.search(r"(?mi)^\s*BC_Type\s*=\s*([^!\r\n]+)", text)
    if not match or [int(value) for value in match.group(1).split()] != [4, 4, 6, 6]:
        fail(f"{path} does not use BC_Type=4 4 6 6")
    check_dis_controls(path)


def check_log(path, expected_steps, max_inner):
    text = path.read_text(errors="replace")
    lowered = text.lower()
    for pattern in ("unsteady inner failure", "step rejected", "rejected", "nonfinite",
                    "mpi_abort", "segmentation fault", "retry"):
        if pattern in lowered:
            fail(f"{path} contains {pattern!r}")
    rows = [(int(step), int(inner)) for step, inner in re.findall(
        r"Unsteady step=(\d+)\s+inner=(\d+)", text
    )]
    expected = list(range(1, expected_steps + 1))
    if [step for step, _ in rows] != expected:
        fail(f"{path} has unexpected accepted steps: {rows[:3]} ... {rows[-3:]}")
    maximum = max(inner for _, inner in rows)
    if maximum >= max_inner:
        fail(f"{path} hit Unsteady_MaxInner={max_inner}: maximum inner={maximum}")
    return maximum, len(rows)


def time_scale_metrics(flow, mesh, config_path, dt, reference_dt):
    nx = config_int(config_path, "Nx_all")
    ny = config_int(config_path, "Ny_all")
    cell_x, cell_y, areas, length1, length2 = mesh
    cells = nx * ny
    gamma = 5.0 / 3.0
    cfl = config_real(config_path, "CFL")
    u_max = config_real(config_path, "U_MAX")
    v_max = config_real(config_path, "V_MAX")
    intervals = []
    upper_bounds = []
    for cell in range(cells):
        rho = flow[cell]
        rho_u = flow[cells + cell]
        rho_v = flow[2 * cells + cell]
        energy = flow[3 * cells + cell]
        u = rho_u / rho
        v = rho_v / rho
        denominator = energy - 0.5 * (rho_u * rho_u + rho_v * rho_v) / rho
        inverse_temperature = 0.5 * rho / ((gamma - 1.0) * denominator)
        sos = math.sqrt(0.5 * gamma / inverse_temperature)
        x_speed = max(u_max, abs(u) + sos)
        y_speed = max(v_max, abs(v) + sos)
        intervals.append(cfl * areas[cell] /
                         (length2[cell] * x_speed + length1[cell] * y_speed))
        upper_bounds.append(cfl * areas[cell] /
                            (length2[cell] * u_max + length1[cell] * v_max))
    dt_s = [min(value, dt) for value in intervals]
    ratios = [value / dt for value in dt_s]
    return {
        "dt_global": dt,
        "dt_s_min": min(dt_s),
        "dt_s_max": max(dt_s),
        "dt_s_over_dt_min": min(ratios),
        "dt_s_over_dt_max": max(ratios),
        "cfl_interval_min": min(intervals),
        "cfl_interval_max": max(intervals),
        "guaranteed_dt_s_upper_bound": max(upper_bounds),
        "guaranteed_dt_s_upper_bound_over_reference_dt": max(upper_bounds) / reference_dt,
    }


def main():
    root = Path(sys.argv[1])
    matrix_path = root / "temporal_matrix.json"
    if not matrix_path.exists():
        fail(f"missing {matrix_path}")
    matrix = json.loads(matrix_path.read_text())
    nx = int(matrix["nx"])
    ny = int(matrix["ny"])
    candidates = [
        ("dt0010", 0.01, 4),
        ("dt0005", 0.005, 8),
        ("dt00025", 0.0025, 16),
        ("dt000125", 0.00125, 32),
    ]
    reference = ("reference", float(matrix["reference_dt"]), int(matrix["reference_steps"]))
    all_specs = [reference] + candidates
    if not math.isclose(reference[1] * reference[2], matrix["final_time"],
                        rel_tol=0.0, abs_tol=1.0e-14):
        fail("reference does not reach final_time exactly")
    for _, dt, steps in candidates:
        if not math.isclose(dt * steps, matrix["final_time"], rel_tol=0.0, abs_tol=1.0e-14):
            fail("candidate does not reach final_time exactly")

    timings = {}
    timing_lines = (root / "timings.tsv").read_text().splitlines()[1:]
    for line in timing_lines:
        fields = line.split("\t")
        if len(fields) != 2:
            fail(f"malformed timing line: {line}")
        timings[fields[0]] = float(fields[1])
    expected_timing_names = {"base", "reference", *(name for name, _, _ in candidates)}
    if set(timings) != expected_timing_names:
        fail(f"incomplete temporal timing records: {sorted(timings)}")
    total_runtime = sum(timings.values())
    if not math.isfinite(total_runtime) or total_runtime > 30.0 * 60.0:
        fail(f"temporal experiment exceeded 30 minutes: {total_runtime:.3f}s", 124)

    smooth_metadata = json.loads((root / "smooth_bump_metadata.json").read_text())
    if not smooth_metadata["moment_compatibility"]["passed"]:
        fail("smooth bump helper did not report moment compatibility")
    initial_flow_bytes = (root / "initial_flow2d.dat").read_bytes()
    initial_vdf_bytes = (root / "initial_vdf2d.dat").read_bytes()
    mesh = read_mesh(root / "initial_Mesh2d.dat", nx, ny)
    initial_flow, initial_vdf, nvel = read_restart(
        root / "initial_flow2d.dat", root / "initial_vdf2d.dat", nx, ny
    )
    compatibility = verify_restart_compatibility(initial_flow, initial_vdf, nx, ny, nvel)
    initial_primitive = primitive_from_flow(initial_flow, nx, ny)

    all_fields = {}
    max_inner = {}
    diagnostic_counts = {}
    dt_s_metrics = {}
    for name, dt, steps in all_specs:
        case_dir = root / name
        if (case_dir / "input_flow2d.dat").read_bytes() != initial_flow_bytes:
            fail(f"{name} did not start from byte-identical smooth flow2d.dat")
        if (case_dir / "input_vdf2d.dat").read_bytes() != initial_vdf_bytes:
            fail(f"{name} did not start from byte-identical smooth vdf2d.dat")
        if (case_dir / "Mesh2d.dat").read_bytes() != (root / "initial_Mesh2d.dat").read_bytes():
            fail(f"{name} did not use the common mesh")
        config_path = case_dir / "Config.ini"
        check_config(config_path, matrix, dt, steps)
        maximum, count = check_log(case_dir / "solver.log", steps, int(matrix["max_inner"]))
        max_inner[name] = maximum
        diagnostic_counts[name] = count
        output = case_dir / f"Flow2d_{steps}.plt"
        if not output.exists():
            fail(f"{name} did not reach final output {output.name}")
        all_fields[name] = read_plt(output, nx, ny)
        dt_s_metrics[name] = time_scale_metrics(
            initial_flow, mesh, config_path, dt, float(matrix["reference_dt"])
        )

    reference_fields = all_fields["reference"]
    for name, values in all_fields.items():
        if any(abs(left - right) > 1.0e-13 for left, right in zip(reference_fields[0], values[0])) or \
                any(abs(left - right) > 1.0e-13 for left, right in zip(reference_fields[1], values[1])):
            fail(f"{name} final output uses a different cell-center grid")

    initial_fields = initial_primitive
    final_perturbations = {}
    for name, values in all_fields.items():
        final_perturbations[name] = [
            [value - initial for value, initial in zip(values[index + 2], initial_fields[index])]
            for index in range(4)
        ]
    reference_perturbation = final_perturbations["reference"]
    field_names = ("rho", "u", "v", "T")
    combined = []
    per_variable = {field_name: {} for field_name in field_names}
    for name, _, _ in candidates:
        candidate_perturbation = final_perturbations[name]
        combined.append(combined_metric(reference_perturbation, candidate_perturbation))
        for index, field_name in enumerate(field_names):
            per_variable[field_name][name] = field_metric(
                reference_perturbation[index], candidate_perturbation[index]
            )

    field_errors = [metric["relative_l2"] for metric in combined]
    field_orders = observed_orders(field_errors)
    numerical_floor = max(
        100.0 * sys.float_info.epsilon * combined_metric(
            reference_perturbation, reference_perturbation
        )["reference_norm"],
        10.0 * float(matrix["inner_tolerance"]),
    )
    monotone = all(left > right for left, right in zip(field_errors, field_errors[1:]))
    last_three_above_floor = all(value > numerical_floor for value in field_errors[-3:])
    last_two_orders_second = len(field_orders) >= 2 and all(
        order >= 1.8 for order in field_orders[-2:]
    )
    formal = monotone and last_three_above_floor and last_two_orders_second

    reference_dt = float(matrix["reference_dt"])
    max_guaranteed_dt_s = max(
        metrics["guaranteed_dt_s_upper_bound"] for metrics in dt_s_metrics.values()
    )
    all_dt_s_less_reference = max_guaranteed_dt_s < reference_dt
    if not all_dt_s_less_reference:
        fail(f"guaranteed dt_s is not below reference_dt: {max_guaranteed_dt_s}")

    metrics = {
        "inputs": matrix,
        "total_runtime_seconds": total_runtime,
        "runtime_seconds": timings,
        "restart_compatibility": compatibility,
        "maximum_inner_iterations": max_inner,
        "implicit_diagnostic_counts": diagnostic_counts,
        "dt_s": dt_s_metrics,
        "all_dt_s_less_reference_dt": all_dt_s_less_reference,
        "reference_perturbation_norm": combined_metric(
            reference_perturbation, reference_perturbation
        )["reference_norm"],
        "combined_field": combined,
        "per_variable": per_variable,
        "field_errors": field_errors,
        "field_orders": field_orders,
        "numerical_floor": numerical_floor,
        "last_three_field_errors_above_floor": last_three_above_floor,
        "field_errors_monotone": monotone,
        "last_two_field_orders_at_least_1_8": last_two_orders_second,
        "formal_second_order_resolved": formal,
    }
    (root / "temporal_order_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print("temporal run root:", root)
    print(json.dumps(metrics, indent=2))
    if not formal:
        raise SystemExit(3)
    print("PASS implicit-unsteady smooth-bump temporal order")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        fail("usage: check_implicit_unsteady_temporal_order.py RUN_ROOT")
    main()
