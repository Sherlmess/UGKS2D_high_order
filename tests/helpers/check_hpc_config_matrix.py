#!/usr/bin/env python3
import math
import re
import sys
from pathlib import Path


ROOT = (
    Path(sys.argv[1]).resolve()
    if len(sys.argv) > 1
    else Path(__file__).resolve().parents[2]
)
CONFIG_DIR = ROOT / "hpc" / "configs"

CASES = {
    "Config.preflight.ini": ("ImplicitSteady", 100.0, 50, 10, 50, 0),
    "Config.explicit.ini": ("Explicit", 100.0, 2_000_000, 100, 20_000, 1),
    "Config.implicit-ncfl-1.ini": ("ImplicitSteady", 1.0, 2_000_000, 100, 20_000, 1),
    "Config.implicit-ncfl-10.ini": ("ImplicitSteady", 10.0, 2_000_000, 100, 20_000, 1),
    "Config.implicit-ncfl-100.ini": ("ImplicitSteady", 100.0, 2_000_000, 100, 20_000, 1),
    "Config.implicit-ncfl-1000.ini": ("ImplicitSteady", 1000.0, 2_000_000, 100, 20_000, 1),
}

COMMON_NUMERIC = {
    "Iflag_InitFlow": 0,
    "Ma": 0.5,
    "Re": 200.0,
    "AOA": 0.0,
    "Tinf": 1.0,
    "Tw": 2.0,
    "ck": 3,
    "Pr": 0.72,
    "flag_vislaw": 0,
    "TimeStep_method": 2,
    "CFL": 0.5,
    "EPS": 1.0e-7,
    "Quadrature_method": 2,
    "Iflag_WriteSurf": 1,
    "Nx_all": 520,
    "Ny_all": 350,
    "Nu": 16,
    "Nv": 16,
    "U_MIN": -3.686007162724397,
    "U_MAX": 3.686007162724397,
    "V_MIN": -3.686007162724397,
    "V_MAX": 3.686007162724397,
    "NP_X": 20,
    "NP_Y": 14,
    "Implicit_ResidualTol": 1.0e-7,
    "Implicit_MomentTol": 1.0e-12,
}


def assignment(text: str, key: str) -> str:
    matches = re.findall(rf"(?mi)^\s*{re.escape(key)}\s*=\s*([^!\r\n]+)", text)
    if len(matches) != 1:
        raise AssertionError(f"{key}: expected one assignment, found {len(matches)}")
    return matches[0].strip()


def number(text: str, key: str) -> float:
    return float(assignment(text, key).replace("d", "e").replace("D", "E"))


def dis_flags(text: str) -> tuple[int, int]:
    blocks = text.split("$DIS_para")
    if len(blocks) < 3:
        raise AssertionError("missing $DIS_para block")
    values = []
    for line in blocks[1].splitlines():
        token = line.split("!", 1)[0].strip()
        if token and not token.startswith("#"):
            values.append(token)
        if len(values) == 2:
            break
    if len(values) != 2:
        raise AssertionError("missing DIS/fringe flags")
    return int(values[0]), int(values[1])


for filename, (method, ncfl, max_iter, show, write, write_restart) in CASES.items():
    path = CONFIG_DIR / filename
    if not path.is_file():
        raise SystemExit(f"missing HPC configuration: {path}")
    text = path.read_text()
    if re.search(r"(?mi)^\s*Reconstruction_method\s*=", text):
        raise SystemExit(
            f"{filename}: Reconstruction_method is not a member of /config_para/"
        )
    for key, expected in COMMON_NUMERIC.items():
        actual = number(text, key)
        if not math.isclose(actual, float(expected), rel_tol=1.0e-12, abs_tol=1.0e-14):
            raise SystemExit(f"{filename}: {key}={actual}, expected {expected}")
    actual_method = assignment(text, "TimeIntegration_method").strip("'\"")
    if actual_method != method:
        raise SystemExit(f"{filename}: method={actual_method}, expected {method}")
    exact = {
        "Implicit_nCFL": ncfl,
        "MAX_ITER": max_iter,
        "Implicit_MaxIter": max_iter,
        "Nstep_show": show,
        "Nstep_write": write,
        "Iflag_WriteRestart": write_restart,
    }
    for key, expected in exact.items():
        if not math.isclose(number(text, key), float(expected), rel_tol=1.0e-12, abs_tol=1.0e-14):
            raise SystemExit(f"{filename}: {key} does not match {expected}")
    if dis_flags(text) != (0, 0):
        raise SystemExit(f"{filename}: DIS and fringe must both be zero")
    if int(number(text, "NP_X") * number(text, "NP_Y")) != 280:
        raise SystemExit(f"{filename}: MPI partition product must be 280")

print("PASS HPC configuration matrix")
