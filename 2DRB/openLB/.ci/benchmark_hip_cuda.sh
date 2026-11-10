#!/usr/bin/env nix
#!nix --extra-experimental-features nix-command --extra-experimental-features flakes develop ../.#hip-cuda --command bash

# Benchmark the GPU_HIP backend with hipcc-style nvcc targeting NVIDIA
# (HIP_PLATFORM=nvidia). Mirrors benchmark_cuda.sh but routes population storage
# and kernels through HIP, so the HIP-on-CUDA path gets the same runtime coverage
# as native CUDA. Sourced entirely from nixpkgs via the hip-cuda dev shell, so it
# runs on the same NVIDIA has-gpu runner without a system ROCm install.

set -euo pipefail

export USE_EMBEDDED_DEPENDENCIES=ON
# The real CUDA driver (/run/opengl-driver/lib, vs the nixpkgs stub) is put on
# LD_LIBRARY_PATH by the .#hip-cuda dev shell itself — no export needed here.
export CUDA_VISIBLE_DEVICES=$(( $(nvidia-smi --list-gpus | wc -l) > 1 ? 1 : 0 ))

echo "Testing on $(nvidia-smi --query-gpu=name --format=csv,noheader --id=$CUDA_VISIBLE_DEVICES) (GPU_HIP -> nvidia, sm_$CUDA_ARCH)"

make -j4

# Note: the gtest-based benchmark tests (e.g. testPoiseuille2DBGK) are not built
# here. Under nvcc -rdc=true the device link step chokes on nixpkgs gtest's
# registration objects; this is a gtest+nvcc toolchain issue independent of the
# HIP-on-CUDA path. Runtime correctness/performance is exercised by the
# cavity3dBenchmark run below.

pushd examples/laminar/cavity3dBenchmark
make
# cavity3dBenchmark runs no statistics in its timed loop, so its output never
# contains NaN/inf even on divergence -- the only meaningful signal is throughput.
# Guard against a broken build/runtime by asserting the averaged MLUPs stays
# within a band of the value expected for this runner (RTX A5000, sm_86).
# Update EXPECTED_MLUPS if the runner hardware changes.
EXPECTED_MLUPS=3029
TOLERANCE=0.10
for run in {1..5}; do
    # reporter::Timer writes each run's "size, steps, ranks, MLUPs" throughput row to
    # profilingData; take its data row (MLUPs is the last comma-separated field).
    ./cavity3d --RESOLUTION 100 --TIME_STEPS 1000 --VTK_ENABLED 0
    tail -n 1 tmp/profilingData/performance.csv | tee --append raw_performance.csv
done
avg=$(awk -F',' '{sum+=$NF; ++n} END { print sum/n }' < raw_performance.csv)
echo "average_cavity3d_hip_cuda_performance_n100 $avg" > performance.txt
awk -v a="$avg" -v e="$EXPECTED_MLUPS" -v t="$TOLERANCE" 'BEGIN{
  lo=e*(1-t); hi=e*(1+t);
  if (a<lo || a>hi) { printf "FAIL: cavity3d avg %.1f MLUPs outside [%.1f, %.1f] (expected %.0f +/-%d%%)\n", a, lo, hi, e, t*100; exit 1 }
  printf "OK: cavity3d avg %.1f MLUPs within [%.1f, %.1f]\n", a, lo, hi
}'
popd
