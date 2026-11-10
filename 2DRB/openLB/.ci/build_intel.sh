#!/usr/bin/env bash

# Intel oneAPI C++ (icpx) compile coverage in the stock public hpckit image
# (intel/oneapi-hpckit). Non-MPI CPU build, no SYCL offload.

set -euo pipefail

if ! command -v icpx >/dev/null 2>&1 && [ -f /opt/intel/oneapi/setvars.sh ]; then
    # shellcheck disable=SC1091
    source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
fi

icpx --version | head -1

NCORES=${NCORES:-4}

cp config/cpu_intel_ci.mk config.mk
make clean-all

make -j"$NCORES"
make -j"$NCORES" samples FEATURES=CIPIPELINE
