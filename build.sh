#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
make -j"${BUILD_JOBS:-4}" FC="${FC:-mpifort}" bin/ugks_solver
