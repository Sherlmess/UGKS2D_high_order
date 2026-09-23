#!/usr/bin/env python3
import argparse
import json
import math
import struct
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


def read_records(path):
    data = path.read_bytes()
    records = []
    offset = 0
    while offset < len(data):
        if offset + 4 > len(data):
            raise ValueError(f"truncated record marker in {path}")
        length = struct.unpack_from("=i", data, offset)[0]
        if length < 0:
            raise ValueError(f"negative record length in {path}")
        start = offset + 4
        end = start + length
        if end + 4 > len(data):
            raise ValueError(f"truncated record payload in {path}")
        trailer = struct.unpack_from("=i", data, end)[0]
        if trailer != length:
            raise ValueError(f"record marker mismatch in {path}")
        records.append(data[start:end])
        offset = end + 4
    return records


def write_records(path, records):
    output = bytearray()
    for payload in records:
        marker = struct.pack("=i", len(payload))
        output.extend(marker)
        output.extend(payload)
        output.extend(marker)
    path.write_bytes(output)


def read_ints(payload, count, label):
    if len(payload) != 4 * count:
        raise ValueError(f"{label} has {len(payload)} bytes, expected {4 * count}")
    return struct.unpack("=" + "i" * count, payload)


def read_doubles(payload, count, label):
    if len(payload) != 8 * count:
        raise ValueError(f"{label} has {len(payload)} bytes, expected {8 * count}")
    values = array("d")
    values.frombytes(payload)
    if len(values) != count:
        raise ValueError(f"{label} has {len(values)} doubles, expected {count}")
    return list(values)


def pack_doubles(values):
    return array("d", values).tobytes()


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


def moments(h, b, u, v, weight):
    return [
        sum(weight[k] * h[k] for k in range(len(h))),
        sum(weight[k] * u[k] * h[k] for k in range(len(h))),
        sum(weight[k] * v[k] * h[k] for k in range(len(h))),
        0.5 * sum(weight[k] * ((u[k] * u[k] + v[k] * v[k]) * h[k] + b[k])
                  for k in range(len(h))),
    ]


def solve_dense4(matrix, rhs):
    augmented = [list(row) + [value] for row, value in zip(matrix, rhs)]
    for column in range(4):
        pivot = max(range(column, 4), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= 1000.0 * 2.220446049250313e-16:
            raise ValueError("singular four-moment projection matrix")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(4):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [left - factor * right
                              for left, right in zip(augmented[row], augmented[column])]
    return [augmented[row][4] for row in range(4)]


def project_to_moments(h, b, target, u, v, weight):
    scale_h = max(max(abs(value) for value in h), 1.0)
    scale_b = max(max(abs(value) for value in b), 1.0)
    kernel_b = [0.0, 0.0, 0.0, 0.5]
    gram = [[0.0 for _ in range(4)] for _ in range(4)]
    current = moments(h, b, u, v, weight)
    rhs = [expected - actual for expected, actual in zip(target, current)]
    for k in range(len(h)):
        kernel_h = [1.0, u[k], v[k], 0.5 * (u[k] * u[k] + v[k] * v[k])]
        weight_h = max(abs(h[k]), 2.220446049250313e-16 * scale_h)
        weight_b = max(abs(b[k]), 2.220446049250313e-16 * scale_b)
        for row in range(4):
            for column in range(4):
                gram[row][column] += weight[k] * (
                    weight_h * kernel_h[row] * kernel_h[column] +
                    weight_b * kernel_b[row] * kernel_b[column]
                )
    lam = solve_dense4(gram, rhs)
    projected_h = list(h)
    projected_b = list(b)
    for k in range(len(h)):
        kernel_h = [1.0, u[k], v[k], 0.5 * (u[k] * u[k] + v[k] * v[k])]
        weight_h = max(abs(h[k]), 2.220446049250313e-16 * scale_h)
        weight_b = max(abs(b[k]), 2.220446049250313e-16 * scale_b)
        projected_h[k] += weight_h * sum(kernel_h[a] * lam[a] for a in range(4))
        projected_b[k] += 0.5 * weight_b * lam[3]
    return projected_h, projected_b


def read_mesh(path):
    records = read_records(path)
    if len(records) != 2:
        raise ValueError(f"{path} must contain two records")
    nxp, nyp = read_ints(records[0], 2, "mesh dimensions")
    node_count = nxp * nyp
    coordinates = read_doubles(records[1], 2 * node_count, "mesh coordinates")
    x_nodes = coordinates[:node_count]
    nx, ny = nxp - 1, nyp - 1
    cell_x = []
    for j in range(ny):
        for i in range(nx):
            lower_left = j * nxp + i
            lower_right = lower_left + 1
            upper_left = (j + 1) * nxp + i
            upper_right = upper_left + 1
            cell_x.append(0.25 * (x_nodes[lower_left] + x_nodes[lower_right] +
                                  x_nodes[upper_left] + x_nodes[upper_right]))
    return nx, ny, cell_x


def read_restart(flow_path, vdf_path):
    flow_records = read_records(flow_path)
    if len(flow_records) != 3:
        raise ValueError(f"{flow_path} must contain three records")
    if read_ints(flow_records[0], 1, "flow header") != (1,):
        raise ValueError(f"unexpected flow header in {flow_path}")
    nx, ny, nvar = read_ints(flow_records[1], 3, "flow dimensions")
    if nvar != 4:
        raise ValueError(f"{flow_path} is not a four-variable restart")
    cells = nx * ny
    flow = read_doubles(flow_records[2], cells * nvar, "flow data")

    vdf_records = read_records(vdf_path)
    if len(vdf_records) != 3:
        raise ValueError(f"{vdf_path} must contain three records")
    if read_ints(vdf_records[0], 1, "vdf header") != (1,):
        raise ValueError(f"unexpected vdf header in {vdf_path}")
    vdf_nx, vdf_ny, vdf_nvar = read_ints(vdf_records[1], 3, "vdf dimensions")
    if (vdf_nx, vdf_ny) != (nx, ny) or vdf_nvar % 2 != 0:
        raise ValueError(f"inconsistent vdf dimensions in {vdf_path}")
    nvel = vdf_nvar // 2
    vdf = read_doubles(vdf_records[2], cells * vdf_nvar, "vdf data")
    if nvel != 16 * 16:
        raise ValueError(f"expected GH16 VDF with 256 velocities, found {nvel}")
    return flow_records, vdf_records, nx, ny, flow, vdf, nvel


def scale_restart(flow_records, vdf_records, nx, ny, flow, vdf, nvel, cell_x,
                  flow_out, vdf_out, metadata_out):
    cells = nx * ny
    u, v, weight = gh16_velocity_space()
    factors = []
    for x in cell_x:
        if 0.25 <= x <= 0.75:
            factor = 1.0 + 1.0e-5 * math.sin(math.pi * (x - 0.25) / 0.5) ** 4
        else:
            factor = 1.0
        factors.append(factor)

    compatible_h = list(vdf[:nvel * cells])
    compatible_b = list(vdf[nvel * cells:])
    base_max_abs_error = 0.0
    base_max_relative_error = 0.0
    for cell in range(cells):
        target = [flow[component * cells + cell] for component in range(4)]
        h = [compatible_h[velocity * cells + cell] for velocity in range(nvel)]
        b = [compatible_b[velocity * cells + cell] for velocity in range(nvel)]
        current = moments(h, b, u, v, weight)
        scale = max(1.0, max(abs(value) for value in target))
        for actual, expected in zip(current, target):
            error = abs(actual - expected)
            base_max_abs_error = max(base_max_abs_error, error)
            base_max_relative_error = max(base_max_relative_error, error / scale)
        projected_h, projected_b = project_to_moments(h, b, target, u, v, weight)
        for velocity in range(nvel):
            compatible_h[velocity * cells + cell] = projected_h[velocity]
            compatible_b[velocity * cells + cell] = projected_b[velocity]

    scaled_flow = list(flow)
    scaled_vdf = compatible_h + compatible_b
    for cell, factor in enumerate(factors):
        for component in range(4):
            scaled_flow[component * cells + cell] *= factor
        for velocity in range(nvel):
            scaled_vdf[velocity * cells + cell] *= factor
            scaled_vdf[(nvel + velocity) * cells + cell] *= factor

    max_abs_error = 0.0
    max_relative_error = 0.0
    for cell in range(cells):
        target = [scaled_flow[component * cells + cell] for component in range(4)]
        h = [scaled_vdf[velocity * cells + cell] for velocity in range(nvel)]
        b = [scaled_vdf[(nvel + velocity) * cells + cell] for velocity in range(nvel)]
        computed = moments(h, b, u, v, weight)
        scale = max(1.0, max(abs(value) for value in target))
        for actual, expected in zip(computed, target):
            error = abs(actual - expected)
            max_abs_error = max(max_abs_error, error)
            max_relative_error = max(max_relative_error, error / scale)

    tolerance = 1.0e-10
    if max_relative_error > tolerance:
        raise ValueError(
            f"smooth bump restart is not discretely compatible: "
            f"max_relative_error={max_relative_error:.6e} > {tolerance:.6e}"
        )

    write_records(flow_out, [flow_records[0], flow_records[1],
                             pack_doubles(scaled_flow)])
    write_records(vdf_out, [vdf_records[0], vdf_records[1],
                            pack_doubles(scaled_vdf)])
    metadata = {
        "nx": nx,
        "ny": ny,
        "nvel": nvel,
        "interval": [0.25, 0.75],
        "amplitude": 1.0e-5,
        "factor_min": min(factors),
        "factor_max": max(factors),
        "modified_cells": sum(factor != 1.0 for factor in factors),
        "production_base_moment_compatibility": {
            "max_absolute_error": base_max_abs_error,
            "max_relative_error": base_max_relative_error,
            "projected_before_smooth_scaling": True,
        },
        "moment_compatibility": {
            "max_absolute_error": max_abs_error,
            "max_relative_error": max_relative_error,
            "tolerance": tolerance,
            "passed": True,
        },
    }
    metadata_out.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main():
    parser = argparse.ArgumentParser(description="make a compatible smooth-bump restart")
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--flow-in", type=Path, required=True)
    parser.add_argument("--vdf-in", type=Path, required=True)
    parser.add_argument("--flow-out", type=Path, required=True)
    parser.add_argument("--vdf-out", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args()

    nx, ny, cell_x = read_mesh(args.mesh)
    flow_records, vdf_records, flow_nx, flow_ny, flow, vdf, nvel = read_restart(
        args.flow_in, args.vdf_in
    )
    if (flow_nx, flow_ny) != (nx, ny):
        raise SystemExit(
            f"mesh/restart dimensions differ: mesh={nx}x{ny}, restart={flow_nx}x{flow_ny}"
        )
    metadata = scale_restart(
        flow_records, vdf_records, nx, ny, flow, vdf, nvel, cell_x,
        args.flow_out, args.vdf_out, args.metadata
    )
    print(json.dumps(metadata, sort_keys=True))


if __name__ == "__main__":
    main()
