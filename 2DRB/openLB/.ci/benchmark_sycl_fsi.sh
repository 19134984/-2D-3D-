#!/usr/bin/env bash

# SYCL FSI benchmarks on the OpenCL CPU device: rigidValve2d angle(t) regression
# + wind-farm turbine-torque regression. Runs lean (ADIOS slices and VTK off).
# Runs in the prebuilt oneAPI image (.ci/docker/sycl-ci.Dockerfile).

set -euo pipefail

# ADIOS2/OpenMPI ship in the prebuilt image; install on demand if missing.
if [ ! -f /usr/include/adios2.h ] || ! command -v /usr/bin/mpicxx >/dev/null 2>&1; then
    apt-get update
    apt-get install -y --no-install-recommends \
        openmpi-bin libopenmpi-dev \
        libadios2-mpi-c++-dev libadios2-mpi-c++11-dev libadios2-mpi-core-dev \
        libadios2-common-c++11-dev libadios2-common-core-dev
fi

# numpy/matplotlib back the valve-profile comparison (compare_valve_profile.py).
if ! python3 -c 'import numpy, matplotlib' >/dev/null 2>&1; then
    apt-get update
    apt-get install -y --no-install-recommends python3-numpy python3-matplotlib
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

ROOT=$PWD

# rigidValve2d: angle(t) regression vs the git-tracked reference profile via
# compare_valve_profile.py (phase-robust functionals + RMS backstop); overlay
# plot published as a CI artifact. Last argument is the RMS backstop (rad).
pushd examples/fsi/rigidValve2d
make -j"$NCORES"
LD_LIBRARY_PATH="$PWD:${LD_LIBRARY_PATH:-}" ./rigidValve2d --RESOLUTION 40 --VTK_ENABLED 0
export MPLCONFIGDIR=$PWD/.mplconfig
python3 "$ROOT/.ci/compare_valve_profile.py" \
    "$ROOT/.ci/reference/rigidValve2d_n40_angle.csv" \
    tmp/gnuplotData/data/valve_n40.csv \
    'gcc-sycl' 3e-2 valve_profile_n40.png
popd

# Wind farm: OLB_MIXED_MODE_INCLUDE_CPP folds wind.cpp into the SYCL device TU.
# Only present in the development repository.
if [ ! -d apps/adrian/fsi/windFarm3d ]; then
    exit 0
fi

pushd apps/adrian/fsi/windFarm3d
make -j"$NCORES" OLB_MIXED_MODE_INCLUDE_CPP=ON

# Lean run (--without-adios-slices; VTK off by default).
LD_LIBRARY_PATH="$PWD:${LD_LIBRARY_PATH:-}" mpirun -n 1 ./wind --resolution 40 --maxPhysT 2 --without-adios-slices

# Turbine-torque regression vs a pinned reference (band absorbs compiler drift).
# Reference reflects the STLreader::getY1 two-voxel wall-distance model.
CSV=tmp/gnuplotData/data/turbine_0.csv
TORQUE=$(tail -n 1 "$CSV" | cut -d, -f3)
echo "windFarm turbine torque = ${TORQUE} (reference -24.7875)"
awk -v v="$TORQUE" -v r=-24.7875 'BEGIN {
  if (v+0 == 0)            { print "FAIL: turbine torque is zero/NaN"; exit 1 }
  d = (v - r) / r; if (d < 0) d = -d;
  if (d > 0.05)            { printf "FAIL: torque %.5f deviates %.2f%% from reference %.5f\n", v, d*100, r; exit 1 }
  printf "OK: turbine torque %.5f within %.2f%% of reference %.5f\n", v, d*100, r
}'
popd
