#!/usr/bin/env bash

# SYCL (icpx) build-only coverage; the runs live in benchmark_sycl. Runs in the
# prebuilt oneAPI image (.ci/docker/sycl-ci.Dockerfile).

set -euo pipefail

# Source setvars in case a slimmer image doesn't put oneAPI on PATH.
if ! command -v icpx >/dev/null 2>&1 && [ -f /opt/intel/oneapi/setvars.sh ]; then
    # shellcheck disable=SC1091
    source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
fi

icpx --version | head -1

NCORES=${NCORES:-4}

cp config/cpu_sycl_ci.mk config.mk
make clean-all
make -j"$NCORES"
make -j"$NCORES" \
     examples/laminar/cavity3dBenchmark/ \
     examples/turbulence/tgv3d/
