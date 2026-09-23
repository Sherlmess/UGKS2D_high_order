#!/usr/bin/env python3
import argparse
import math
import re
from pathlib import Path


def fields(path: Path) -> list[list[float]]:
    lines = path.read_text().replace("D", "E").splitlines()
    zone_index = next(i for i, line in enumerate(lines) if "ZONE" in line.upper())
    match_i = re.search(r"I\s*=\s*(\d+)", lines[zone_index], re.I)
    match_j = re.search(r"J\s*=\s*(\d+)", lines[zone_index], re.I)
    if match_i is None or match_j is None:
        raise SystemExit(f"cannot read zone dimensions: {path}")
    count = int(match_i.group(1)) * int(match_j.group(1))
    values = [float(token) for line in lines[zone_index + 1:] for token in line.split()]
    if len(values) < 12 * count or not all(math.isfinite(value) for value in values):
        raise SystemExit(f"invalid block data: {path}")
    blocks = [values[index * count:(index + 1) * count] for index in range(12)]
    return blocks[2:6]


def l2_relative(candidate: list[float], reference: list[float]) -> float:
    numerator = math.sqrt(sum((a - b) ** 2 for a, b in zip(candidate, reference)))
    denominator = max(math.sqrt(sum(value * value for value in reference)), 1.0e-300)
    return numerator / denominator


def linf_relative(candidate: list[float], reference: list[float]) -> float:
    numerator = max(abs(a - b) for a, b in zip(candidate, reference))
    denominator = max(max(abs(value) for value in reference), 1.0e-300)
    return numerator / denominator


parser = argparse.ArgumentParser()
parser.add_argument("reference", type=Path)
parser.add_argument("candidate", type=Path)
parser.add_argument("--reference-surface", type=Path)
parser.add_argument("--candidate-surface", type=Path)
args = parser.parse_args()
reference_fields = fields(args.reference)
candidate_fields = fields(args.candidate)
names = ["rho", "u", "v", "temperature"]
for name, reference, candidate in zip(names, reference_fields, candidate_fields):
    if len(reference) != len(candidate):
        raise SystemExit(f"field dimensions differ for {name}: {len(reference)} != {len(candidate)}")
    reference_scale = max(abs(value) for value in reference)
    absolute_linf = max(abs(a - b) for a, b in zip(candidate, reference))
    if reference_scale < 1.0e-8:
        print(f"{name}: absLinf={absolute_linf:.8e}")
        if absolute_linf >= 1.0e-8:
            raise SystemExit(f"near-zero field comparison failed for {name}")
    else:
        error_l2 = l2_relative(candidate, reference)
        error_linf = linf_relative(candidate, reference)
        print(f"{name}: L2={error_l2:.8e} Linf={error_linf:.8e}")
        if error_l2 >= 1.0e-6 or error_linf >= 1.0e-4:
            raise SystemExit(f"field comparison failed for {name}")
if (args.reference_surface is None) != (args.candidate_surface is None):
    raise SystemExit("both surface files must be supplied together")
if args.reference_surface is not None:
    def surface_columns(path: Path) -> list[list[float]]:
        rows = [[float(token) for token in line.replace("D", "E").split()]
                for line in path.read_text().splitlines()
                if line.strip() and not line.lstrip().lower().startswith("variables")]
        if not rows or not all(len(row) == 5 and all(math.isfinite(value) for value in row) for row in rows):
            raise SystemExit(f"invalid surface data: {path}")
        return [[row[column] for row in rows] for column in (2, 3, 4)]

    for name, reference, candidate in zip(
            ["Cp", "Cf", "Ch"], surface_columns(args.reference_surface), surface_columns(args.candidate_surface)):
        if len(reference) != len(candidate):
            raise SystemExit(f"surface dimensions differ for {name}: {len(reference)} != {len(candidate)}")
        reference_scale = max(abs(value) for value in reference)
        absolute_linf = max(abs(a - b) for a, b in zip(candidate, reference))
        if reference_scale < 1.0e-5:
            print(f"{name}: absLinf={absolute_linf:.8e}")
            if absolute_linf >= 1.0e-5:
                raise SystemExit(f"near-zero surface comparison failed for {name}")
        else:
            error_l2 = l2_relative(candidate, reference)
            error_linf = linf_relative(candidate, reference)
            print(f"{name}: L2={error_l2:.8e} Linf={error_linf:.8e} absLinf={absolute_linf:.8e}")
            if error_l2 >= 1.0e-6 or error_linf >= 1.0e-4:
                raise SystemExit(f"surface comparison failed for {name}")
print("PASS steady field comparison")
