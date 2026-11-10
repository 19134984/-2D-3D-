#!/usr/bin/env bash

# SYCL smoke test on the OpenCL CPU device: builds the library and runs cavity3d
# (crash/divergence fails via the exit code). FSI regressions (rigidValve2d, wind
# farm) live in benchmark_sycl_fsi.sh. Runs in the prebuilt oneAPI image
# (.ci/docker/sycl-ci.Dockerfile).

set -euo pipefail

# ADIOS2/OpenMPI ship in the prebuilt image; install on demand if missing.
if [ ! -f /usr/include/adios2.h ] || ! command -v /usr/bin/mpicxx >/dev/null 2>&1; then
    apt-get update
    apt-get install -y --no-install-recommends \
        openmpi-bin libopenmpi-dev \
        libadios2-mpi-c++-dev libadios2-mpi-c++11-dev libadios2-mpi-core-dev \
        libadios2-common-c++11-dev libadios2-common-core-dev
fi

# Use the OpenMPI wrappers (matching ADIOS2's MPI) wrapping the Intel compilers.
export PATH=/usr/bin:$PATH
export OMPI_CXX=icpx
export OMPI_CC=icx
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

if ! command -v icpx >/dev/null 2>&1 && [ -f /opt/intel/oneapi/setvars.sh ]; then
    # shellcheck disable=SC1091
    source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
fi

echo "mpicxx => $(mpicxx --showme:version 2>&1 | head -1)"
icpx --version | head -1

NCORES=${NCORES:-4}
export ONEAPI_DEVICE_SELECTOR=opencl:cpu

cp config/cpu_sycl_adios_ci.mk config.mk
make clean-all
make -j"$NCORES"

# cavity3d smoke test (crash/divergence fails via the exit code).
make -j"$NCORES" examples/laminar/cavity3dBenchmark/
pushd examples/laminar/cavity3dBenchmark
LD_LIBRARY_PATH="$PWD:${LD_LIBRARY_PATH:-}" ./cavity3d --RESOLUTION 30 --TIME_STEPS 200 --VTK_ENABLED 0
popd
