#!/usr/bin/env bash

# Intel oneAPI C++ (icpx) with MPI compile coverage in the stock public hpckit
# image (intel/oneapi-hpckit), using Intel MPI's icpx wrapper mpiicpx.

set -euo pipefail

if ! command -v icpx >/dev/null 2>&1 && [ -f /opt/intel/oneapi/setvars.sh ]; then
    # shellcheck disable=SC1091
    source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
fi

mpiicpx --version | head -1

NCORES=${NCORES:-4}

cp config/cpu_intel_mpi_ci.mk config.mk
make clean-all

make -j"$NCORES"
make -j"$NCORES" samples FEATURES=CIPIPELINE
