#!/usr/bin/env bash
# Builds the communication integrity oracle (tests/communication) and runs it
# across MPI rank counts and the directional canary. The toolchain is selected
# by the caller's nix shell. Set COMM_TEST_BUILD_ONLY=1 for compile-only
# coverage where no runnable device is available (HIP on AMD runners without a
# GPU). Set COMM_TEST_DIRECT=1 to run the binary directly instead of via mpirun
# (non-MPI configs whose image nonetheless ships an mpirun, e.g. GPU_SYCL on the
# OpenCL CPU runtime with ONEAPI_DEVICE_SELECTOR=opencl:cpu).

set -uo pipefail
cd "$(dirname "$0")/.."

if command -v nvidia-smi >/dev/null 2>&1; then
  export USE_EMBEDDED_DEPENDENCIES=ON
  # $CUDA_PATH/lib carries libcudart for the HIP-on-CUDA (HIP_PLATFORM=nvidia)
  # single-compiler link; /run/opengl-driver/lib carries the real driver.
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:${CUDA_PATH:+$CUDA_PATH/lib:}/run/opengl-driver/lib
  export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n1 | grep -o '[0-9]' | tr -d '\n')
  export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}
fi

pushd tests/communication
make

if [[ "${COMM_TEST_BUILD_ONLY:-0}" == "1" ]]; then
  echo "COMM_TEST_BUILD_ONLY=1: compile-only, skipping run"
  popd
  exit 0
fi

status=0

if [[ "${COMM_TEST_DIRECT:-0}" != "1" ]] && command -v mpirun >/dev/null 2>&1; then
  run() {
    local ranks=$1; shift
    echo "== mpirun -n ${ranks} $* =="
    LD_LIBRARY_PATH="${PWD}:${LD_LIBRARY_PATH:-}" \
      mpirun -n "${ranks}" -x LD_LIBRARY_PATH ${ONEAPI_DEVICE_SELECTOR:+-x ONEAPI_DEVICE_SELECTOR} "$@" ./testCommIntegrity || status=1
  }
  run 1
  run 2
  run 1 -x OLB_DIRECTIONAL_CANARY=1
  run 2 -x OLB_DIRECTIONAL_CANARY=1
else
  run() {
    echo "== $* ./testCommIntegrity =="
    LD_LIBRARY_PATH="${PWD}:${LD_LIBRARY_PATH:-}" env "$@" ./testCommIntegrity || status=1
  }
  run
  run OLB_DIRECTIONAL_CANARY=1
fi

popd
exit ${status}
