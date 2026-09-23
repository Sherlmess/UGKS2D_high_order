#!/usr/bin/env python3
import subprocess
from pathlib import Path


root = Path(__file__).resolve().parents[2]
text = (root / "makefile").read_text()
required = [
    "ifeq ($(FC), mpiifort)",
    "MOD_FLAG = -module $(OBJ_DIR)",
    "DEBUG_FLAGS = -g -O0 -check all -traceback -fpe0",
    "RELEASE_FLAGS = -O3 -xHost -ipo",
]

for item in required:
    if item not in text:
        raise SystemExit(f"missing mpiifort Makefile contract: {item}")

branch = text.split("ifeq ($(FC), mpiifort)", 1)[1].split("endif", 1)[0]
if "-qopenmp" in branch:
    raise SystemExit("mpiifort release must remain pure MPI")

result = subprocess.run(
    ["make", "-n", "FC=mpiifort", "release"],
    cwd=root,
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=True,
)
if "mpiifort" not in result.stdout or "-O3 -xHost -ipo" not in result.stdout:
    raise SystemExit("mpiifort release dry-run does not expose Intel release flags")

print("PASS mpiifort Makefile contract")
