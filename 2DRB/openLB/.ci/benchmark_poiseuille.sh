#!/usr/bin/env bash
# Builds core and runs the poiseuille physics/serialization benchmarks
# (tests/benchmarks) for the toolchain of the invoking environment. Reused
# across platforms to give directional-streaming physics + save/load coverage
# beyond the communication integrity oracle.
set -euo pipefail

cd "$(dirname "$0")/.."

# On a CUDA runner (nvcc sole-driver shell) select the device and its arch.
if command -v nvidia-smi >/dev/null 2>&1; then
  export USE_EMBEDDED_DEPENDENCIES=ON
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:/run/opengl-driver/lib
  export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n1 | grep -o '[0-9]' | tr -d '\n')
  export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0}
fi

make -j"${NCORES:-8}"

# `make test<X>` resolves to `run-test<X>` in tests/benchmarks (build + run).
# Run sequentially so the serializable cases do not race on their checkpoint files.
for bench in testPoiseuille2DBGK testPoiseuille3DBGK \
             testSerializablePoiseuille2DBGK testSerializablePoiseuille3DBGK; do
  make "${bench}"
done
