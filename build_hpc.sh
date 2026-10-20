#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
module purge
module add compiler/intel/2021.3.0 mpi/intelmpi/2021.3.0
make -j"${BUILD_JOBS:-4}" FC=mpiifort COMPILER=intel CC=icc bin/ugks_solver
