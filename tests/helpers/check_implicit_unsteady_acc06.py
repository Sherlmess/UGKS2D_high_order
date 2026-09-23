#!/usr/bin/env python3
import json
import math
import re
import struct
import sys
from pathlib import Path


FLOAT_RE = r"[+\-0-9.eEdD]+"


def as_float(token):
    return float(token.replace("D", "E").replace("d", "e"))


def read_fortran_records(path):
    raw = path.read_bytes()
    records = []
    offset = 0
    while offset < len(raw):
        if offset + 4 > len(raw):
            raise SystemExit(f"truncated record marker in {path}")
        length = struct.unpack_from("<i", raw, offset)[0]
        offset += 4
        if length < 0 or offset + length + 4 > len(raw):
            raise SystemExit(f"invalid four-byte record marker in {path}: {length}")
        payload = raw[offset : offset + length]
        offset += length
        trailing = struct.unpack_from("<i", raw, offset)[0]
        offset += 4
        if trailing != length:
            raise SystemExit(f"mismatched record markers in {path}: {length} != {trailing}")
        records.append(payload)
    return records


def unpack_ints(payload, count, path):
    if len(payload) != 4 * count:
        raise SystemExit(f"unexpected integer record size in {path}: {len(payload)}")
    return struct.unpack("<" + "i" * count, payload)


def unpack_reals(payload, count, path):
    if len(payload) != 8 * count:
        raise SystemExit(f"unexpected real record size in {path}: {len(payload)}")
    return list(struct.unpack("<" + "d" * count, payload))


def read_restart(path, nx, ny, nvel, kind):
    records = read_fortran_records(path)
    if len(records) != 3:
        raise SystemExit(f"expected three records in {path}, got {len(records)}")
    if unpack_ints(records[0], 1, path) != (1,):
        raise SystemExit(f"unexpected restart header in {path}")
    ncell = nx * ny
    if kind == "flow":
        dims = unpack_ints(records[1], 3, path)
        if dims != (nx, ny, 4):
            raise SystemExit(f"unexpected flow dimensions in {path}: {dims}")
        flat = unpack_reals(records[2], 4 * ncell, path)
        return [flat[k * ncell : (k + 1) * ncell]
                for k in range(4)]
    dims = unpack_ints(records[1], 3, path)
    if dims != (nx, ny, 2 * nvel):
        raise SystemExit(f"unexpected VDF dimensions in {path}: {dims}")
    flat = unpack_reals(records[2], 2 * ncell * nvel, path)
    block = ncell * nvel
    return [flat[:block], flat[block:]]


def numeric_tokens(path):
    values = []
    for token in path.read_text().split():
        try:
            values.append(as_float(token.rstrip(",")))
        except ValueError:
            pass
    return values


def latest_flow_file(directory):
    files = list(directory.glob("Flow2d_*.plt"))
    if not files:
        raise SystemExit(f"missing Flow2d_N.plt in {directory}")
    numbered = []
    for path in files:
        match = re.search(r"_(\d+)\.plt$", path.name)
        if match:
            numbered.append((int(match.group(1)), path))
    if not numbered:
        raise SystemExit(f"missing numbered Flow2d_N.plt in {directory}")
    return max(numbered)[1]


def read_flow_plt(path, nx, ny):
    ncell = nx * ny
    values = numeric_tokens(path)
    expected = 12 * ncell
    if len(values) < expected:
        raise SystemExit(f"not enough numeric data in {path}: {len(values)} < {expected}")
    data = values[-expected:]
    return [data[k * ncell : (k + 1) * ncell] for k in range(12)]


def read_surface(path):
    if not path.exists():
        raise SystemExit(f"missing wall output: {path}")
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) != 5:
            continue
        try:
            row = [as_float(value) for value in fields]
        except ValueError:
            continue
        rows.append(row)
    if not rows:
        raise SystemExit(f"missing numeric wall output rows: {path}")
    columns = list(zip(*rows))
    return [list(column) for column in columns]


def norm_metrics(reference, candidate):
    if len(reference) != len(candidate):
        raise SystemExit("field lengths differ between explicit and implicit outputs")
    if not all(math.isfinite(value) for value in reference + candidate):
        raise SystemExit("nonfinite field in explicit/implicit comparison")
    diff2 = sum((value - ref) ** 2 for ref, value in zip(reference, candidate))
    ref2 = sum(ref ** 2 for ref in reference)
    ref_norm = math.sqrt(ref2)
    relative = math.sqrt(diff2) / max(ref_norm, 1.0e-300)
    absolute = max(abs(value - ref) for ref, value in zip(reference, candidate))
    return {
        "relative_l2": relative,
        "absolute_linf": absolute,
        "reference_norm": ref_norm,
        "near_zero_reference": ref_norm <= 1.0e-14,
    }


def check_gate(label, metric, relative_limit, absolute_limit, either=False):
    relative_ok = metric["relative_l2"] <= relative_limit
    absolute_ok = metric["absolute_linf"] <= absolute_limit
    if metric["near_zero_reference"]:
        relative_ok = True
    passed = (relative_ok or absolute_ok) if either else (relative_ok and absolute_ok)
    if not passed:
        raise SystemExit(
            f"{label} exceeds gate: relative={metric['relative_l2']:.6e}, "
            f"absolute={metric['absolute_linf']:.6e}"
        )


def config_value(path, key):
    pattern = re.compile(rf"(?mi)^\s*{re.escape(key)}\s*=\s*([^!\r\n]+)")
    match = pattern.search(path.read_text())
    if not match:
        raise SystemExit(f"missing {key} in {path}")
    return match.group(1).strip().split()[0].rstrip(",/")


def config_real(path, key):
    return as_float(config_value(path, key))


def config_int(path, key):
    return int(config_value(path, key))


def parse_mesh(path, nx, ny):
    records = read_fortran_records(path)
    if len(records) != 2:
        raise SystemExit(f"unexpected mesh record count in {path}")
    nxp, nyp = unpack_ints(records[0], 2, path)
    if (nxp, nyp) != (nx + 1, ny + 1):
        raise SystemExit(f"unexpected mesh dimensions: {(nxp, nyp)}")
    coords = unpack_reals(records[1], 2 * nxp * nyp, path)
    npoint = nxp * nyp
    return coords[:npoint], coords[npoint:], nxp


def cell_geometry(xs, ys, nxp, nx, ny, i, j):
    p00 = j * nxp + i
    p10 = p00 + 1
    p01 = p00 + nxp
    p11 = p01 + 1
    len_x0 = math.hypot(xs[p10] - xs[p00], ys[p10] - ys[p00])
    len_x1 = math.hypot(xs[p11] - xs[p01], ys[p11] - ys[p01])
    len_y0 = math.hypot(xs[p01] - xs[p00], ys[p01] - ys[p00])
    len_y1 = math.hypot(xs[p11] - xs[p10], ys[p11] - ys[p10])
    area = abs((xs[p00] - xs[p11]) * (ys[p10] - ys[p01]) -
               (xs[p10] - xs[p01]) * (ys[p00] - ys[p11])) * 0.5
    return 0.5 * (len_x0 + len_x1), 0.5 * (len_y0 + len_y1), area


def primitive_from_conserved(w, gamma):
    rho, mom_u, mom_v, energy = w
    velocity_u = mom_u / rho
    velocity_v = mom_v / rho
    denom = energy - 0.5 * (mom_u * mom_u + mom_v * mom_v) / rho
    temperature = 2.0 * (gamma - 1.0) * denom / rho
    return rho, velocity_u, velocity_v, temperature


def cfl_interval_min(flow, xs, ys, nxp, nx, ny, cfl, umax, vmax, gamma):
    ncell = nx * ny
    minimum = math.inf
    for j in range(ny):
        for i in range(nx):
            cell = j * nx + i
            rho, velocity_u, velocity_v, temperature = primitive_from_conserved(
                [flow[k][cell] for k in range(4)], gamma
            )
            if not all(math.isfinite(value) and value > 0.0 for value in (rho, temperature)):
                raise SystemExit("nonphysical state while reconstructing time-scale evidence")
            sound_speed = math.sqrt(0.5 * gamma / temperature)
            wave_u = max(umax, abs(velocity_u) + sound_speed)
            wave_v = max(vmax, abs(velocity_v) + sound_speed)
            length_x, length_y, area = cell_geometry(xs, ys, nxp, nx, ny, i, j)
            interval = cfl * area / (length_y * wave_u + length_x * wave_v)
            minimum = min(minimum, interval)
    return minimum


def check_dimensions_and_config(explicit_dir, implicit_dir):
    if (explicit_dir / "Mesh2d.dat").read_bytes() != (implicit_dir / "Mesh2d.dat").read_bytes():
        raise SystemExit("explicit and implicit meshes are not byte-identical")
    keys = {
        "Nx_all": 8, "Ny_all": 6, "Nu": 16, "Nv": 16, "NP_X": 1, "NP_Y": 1,
        "Quadrature_method": 2, "Iflag_InitFlow": 0, "Iflag_WriteSurf": 1,
        "Iflag_WriteRestart": 1, "MAX_ITER": 1, "Nstep_show": 1, "Nstep_write": 1,
        "TimeStep_method": 0, "FluxWall_type": 2,
    }
    for directory in (explicit_dir, implicit_dir):
        config = directory / "Config.ini"
        for key, expected in keys.items():
            if config_int(config, key) != expected:
                raise SystemExit(f"{directory.name}: {key} is not {expected}")
        if not math.isclose(config_real(config, "ACC"), 0.6, rel_tol=0.0, abs_tol=1.0e-14):
            raise SystemExit(f"{directory.name}: ACC is not 0.6")
        if not math.isclose(config_real(config, "dt_global"), 1.0e-6, rel_tol=0.0, abs_tol=1.0e-15):
            raise SystemExit(f"{directory.name}: dt_global is not 1.0e-6")
        if not math.isclose(config_real(config, "Unsteady_Theta"), 0.5, rel_tol=0.0, abs_tol=1.0e-14):
            raise SystemExit(f"{directory.name}: Unsteady_Theta is not 0.5")
        if not math.isclose(config_real(config, "Unsteady_InnerTol"), 1.0e-10, rel_tol=0.0, abs_tol=1.0e-15):
            raise SystemExit(f"{directory.name}: Unsteady_InnerTol is not 1.0e-10")
        if config_int(config, "Unsteady_MaxInner") != 80:
            raise SystemExit(f"{directory.name}: Unsteady_MaxInner is not 80")
    if config_value(explicit_dir / "Config.ini", "TimeIntegration_method") != "'Explicit'":
        raise SystemExit("explicit run did not use TimeIntegration_method='Explicit'")
    if config_value(implicit_dir / "Config.ini", "TimeIntegration_method") != "'ImplicitUnsteady'":
        raise SystemExit("implicit run did not use TimeIntegration_method='ImplicitUnsteady'")


def check_log(directory, require_unsteady_step):
    log = (directory / "solver.log").read_text(errors="replace")
    lower = log.lower()
    for marker in ("unsteady inner failure", "rejected", "error"):
        if marker in lower:
            raise SystemExit(f"{directory.name} log contains {marker!r}")
    if not require_unsteady_step:
        return log, []
    rows = re.findall(r"Unsteady step=(\d+) inner=(\d+)", log)
    if len(rows) != 1 or int(rows[0][0]) != 1:
        raise SystemExit(f"expected one accepted unsteady step in implicit log: {rows}")
    if int(rows[0][1]) > 80:
        raise SystemExit(f"implicit inner iterations exceeded 80: {rows}")
    return log, rows


def main(root):
    explicit_dir = root / "explicit"
    implicit_dir = root / "implicit"
    check_dimensions_and_config(explicit_dir, implicit_dir)
    explicit_log, _ = check_log(explicit_dir, False)
    implicit_log, rows = check_log(implicit_dir, True)

    nx, ny, nvel = 8, 6, 16 * 16
    flow_explicit = read_restart(explicit_dir / "flow2d.dat", nx, ny, nvel, "flow")
    flow_implicit = read_restart(implicit_dir / "flow2d.dat", nx, ny, nvel, "flow")
    vdf_explicit = read_restart(explicit_dir / "vdf2d.dat", nx, ny, nvel, "vdf")
    vdf_implicit = read_restart(implicit_dir / "vdf2d.dat", nx, ny, nvel, "vdf")

    final_flow_explicit = read_flow_plt(latest_flow_file(explicit_dir), nx, ny)
    final_flow_implicit = read_flow_plt(latest_flow_file(implicit_dir), nx, ny)
    surface_explicit = read_surface(explicit_dir / "Surf_CpCfCh.plt")
    surface_implicit = read_surface(implicit_dir / "Surf_CpCfCh.plt")
    if len(surface_explicit[0]) != nx or len(surface_implicit[0]) != nx:
        raise SystemExit("wall output does not contain exactly Nx rows")

    metrics = {"flow2d_conserved": {}, "vdf": {}, "flow_plt": {}, "surface": {}}
    names = ("rho", "rhoU", "rhoV", "energy")
    for name, reference, candidate in zip(names, flow_explicit, flow_implicit):
        metric = norm_metrics(reference, candidate)
        check_gate(name, metric, 3.0e-4, 2.0e-8)
        metrics["flow2d_conserved"][name] = metric
    for name, reference, candidate in zip(("h", "b"), vdf_explicit, vdf_implicit):
        metric = norm_metrics(reference, candidate)
        check_gate(name, metric, 3.0e-4, 2.0e-8)
        metrics["vdf"][name] = metric
    for name, reference, candidate in zip(("rho", "u", "v", "T"),
                                          final_flow_explicit[2:6], final_flow_implicit[2:6]):
        metric = norm_metrics(reference, candidate)
        check_gate(name, metric, 3.0e-4, 2.0e-8)
        metrics["flow_plt"][name] = metric
    for name, reference, candidate in zip(("Cp", "Cf", "Ch"),
                                          surface_explicit[2:5], surface_implicit[2:5]):
        metric = norm_metrics(reference, candidate)
        check_gate(name, metric, 3.0e-4, 2.0e-6, either=True)
        metrics["surface"][name] = metric

    implicit_config = implicit_dir / "Config.ini"
    xs, ys, nxp = parse_mesh(implicit_dir / "Mesh2d.dat", nx, ny)
    cfl = config_real(implicit_config, "CFL")
    umax = config_real(implicit_config, "U_MAX")
    vmax = config_real(implicit_config, "V_MAX")
    ma = config_real(implicit_config, "Ma")
    tinf = config_real(implicit_config, "Tinf")
    ck = config_int(implicit_config, "ck")
    gamma = (ck + 4.0) / (ck + 2.0)
    initial_u = ma * math.sqrt(gamma / 2.0)
    initial_flow = [[1.0] * (nx * ny), [initial_u] * (nx * ny),
                    [0.0] * (nx * ny), [0.5 / (gamma - 1.0) + 0.5 * initial_u ** 2] * (nx * ny)]
    for cell in range(nx * ny):
        initial_flow[3][cell] = 0.5 / ((gamma - 1.0) * tinf) + 0.5 * initial_u ** 2
    cfl_initial = cfl_interval_min(initial_flow, xs, ys, nxp, nx, ny, cfl, umax, vmax, gamma)
    cfl_final = cfl_interval_min(flow_implicit, xs, ys, nxp, nx, ny, cfl, umax, vmax, gamma)
    dt_global_match = re.search(r"dt_global\s*=\s*(%s)" % FLOAT_RE, implicit_log)
    if not dt_global_match:
        raise SystemExit("implicit solver log has no printed dt_global evidence")
    dt_global = as_float(dt_global_match.group(1))
    if not cfl_initial > dt_global or not cfl_final > dt_global:
        raise SystemExit(
            f"CappedFluxTime evidence does not prove dt_s=dt_global: "
            f"initial_cfl={cfl_initial:.16e}, final_cfl={cfl_final:.16e}, dt={dt_global:.16e}"
        )
    dt_s = min(cfl_initial, dt_global)
    ratio = dt_s / dt_global
    if not math.isclose(ratio, 1.0, rel_tol=0.0, abs_tol=1.0e-14):
        raise SystemExit(f"unexpected dt_s/dt_global ratio: {ratio}")
    metrics["time_scale_evidence"] = {
        "dt_global_printed": dt_global,
        "initial_cfl_interval_min": cfl_initial,
        "final_cfl_interval_min": cfl_final,
        "dt_s_reconstructed": dt_s,
        "dt_s_over_dt_global": ratio,
        "method": "printed dt_global plus independent CappedFluxTime reconstruction from Mesh2d.dat and flow2d.dat",
    }
    metrics["maximum_inner_iterations"] = int(rows[0][1])
    metrics["limits"] = {
        "field_relative_l2": 3.0e-4,
        "field_absolute_linf": 2.0e-8,
        "surface_relative_l2": 3.0e-4,
        "surface_absolute_linf": 2.0e-6,
    }
    (root / "acc06_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

    print("PASS implicit-unsteady ACC=0.6 explicit reduction")
    for group, group_metrics in metrics.items():
        if isinstance(group_metrics, dict) and group not in ("time_scale_evidence", "limits"):
            for name, metric in group_metrics.items():
                print(f"{group}.{name}: relative L2={metric['relative_l2']:.6e} absolute Linf={metric['absolute_linf']:.6e}")
    print(f"maximum inner iterations: {metrics['maximum_inner_iterations']}")
    print("confirmed dt_s/dt_global evidence:", json.dumps(metrics["time_scale_evidence"], sort_keys=True))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_implicit_unsteady_acc06.py RUN_ROOT")
    main(Path(sys.argv[1]))
