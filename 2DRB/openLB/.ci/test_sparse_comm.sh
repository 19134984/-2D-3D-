#!/usr/bin/env bash
# Builds the sparse-field communication oracle (tests/communication/sparse) and
# runs it across MPI rank counts. The toolchain is selected by the caller's nix
# shell. Set COMM_TEST_BUILD_ONLY=1 for compile-only coverage where no runnable
# device is available. Set COMM_TEST_DIRECT=1 to run the binary directly instead
# of via mpirun.

set -uo pipefail
cd "$(dirname "$0")/.."

if command -v nvidia-smi >/dev/null 2>&1; then
  export USE_EMBEDDED_DEPENDENCIES=ON
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:${CUDA_PATH:+$CUDA_PATH/lib:}/run/opengl-driver/lib
  export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n1 | grep -o '[0-9]' | tr -d '\n')
  export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}
fi

pushd tests/communication/sparse
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
      mpirun -n "${ranks}" -x LD_LIBRARY_PATH ${ONEAPI_DEVICE_SELECTOR:+-x ONEAPI_DEVICE_SELECTOR} "$@" ./testSparseComm || status=1
  }
  run 1
  run 2
else
  run() {
    echo "== $* ./testSparseComm =="
    LD_LIBRARY_PATH="${PWD}:${LD_LIBRARY_PATH:-}" env "$@" ./testSparseComm || status=1
  }
  run
fi

popd
exit ${status}
