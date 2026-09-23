#!/usr/bin/env python3
"""Guard the Intel-MPI-compatible import used by the implicit boundary module."""

from pathlib import Path
import re


source = Path("src/mod_ImplicitBoundary.f90").read_text(encoding="utf-8")

if re.search(r"(?im)^\s*use\s+mpi\s*,\s*only\s*:.*\bMPI_Allreduce\b", source):
    raise SystemExit(
        "mod_ImplicitBoundary must not import MPI_Allreduce through a USE MPI ONLY "
        "list; use the project-wide collective import pattern for Intel MPI compatibility"
    )

if not re.search(r"(?im)^\s*use\s+mpi\s*(?:!.*)?$", source):
    raise SystemExit("mod_ImplicitBoundary must use the unrestricted MPI module")

print("Intel MPI module import compatibility check passed")
