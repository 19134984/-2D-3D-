#!/usr/bin/env nix
#!nix --extra-experimental-features nix-command --extra-experimental-features flakes develop ../.#cuda --command bash

set -euo pipefail

export USE_EMBEDDED_DEPENDENCIES=ON
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
export CUDA_VISIBLE_DEVICES=$(( $(nvidia-smi --list-gpus | wc -l) > 1 ? 1 : 0 ))
export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader --id=$CUDA_VISIBLE_DEVICES | grep -o [0-9] | tr -d '\n')

echo "Testing on $(nvidia-smi --query-gpu=name --format=csv,noheader --id=$CUDA_VISIBLE_DEVICES)"

make -j4
make testPoiseuille2DBGK

# Benchmarks run no statistics in the timed loop, so throughput is the only signal.
# Assert the averaged MLUPs stays within TOLERANCE of the value expected for this
# runner (RTX A5000, sm_86); update the expected values below if the runner changes.
# make clean per build because FEATURES=FP16 is not tracked as a make dependency.
TOLERANCE=0.10
: > performance.txt

# TGV is triply-periodic, so it needs a larger domain than the cavity to amortise
# the per-step periodic halo exchange and reach its saturated (bandwidth-bound)
# throughput -- where, having no boundary conditions, it out-runs the cavity.
bench() { # name dir bin make-flags resolution steps expected-MLUPs
  local name=$1 dir=$2 bin=$3 flags=$4 res=$5 steps=$6 expected=$7 csv
  ( cd "$dir" && make clean >/dev/null && make FLOATING_POINT_TYPE=float $flags -j8 )
  csv="$dir/raw_${name}.csv"; : > "$csv"
  for run in {1..5}; do
    ( cd "$dir" && ./"$bin" --RESOLUTION "$res" --TIME_STEPS "$steps" --VTK_ENABLED 0 )
    # reporter::Timer writes each run's throughput row to profilingData; MLUPs is its last field
    tail -n 1 "$dir/tmp/profilingData/performance.csv" | tee --append "$csv"
  done
  local avg; avg=$(awk -F',' '{sum+=$NF; ++n} END { print sum/n }' "$csv")
  echo "average_${name} $avg" >> performance.txt
  awk -v a="$avg" -v e="$expected" -v t="$TOLERANCE" -v L="$name" 'BEGIN{
    lo=e*(1-t); hi=e*(1+t);
    if (a<lo || a>hi) { printf "FAIL: %s avg %.1f MLUPs outside [%.1f, %.1f] (expected %.0f +/-%d%%)\n", L, a, lo, hi, e, t*100; exit 1 }
    printf "OK: %s avg %.1f MLUPs within [%.1f, %.1f]\n", L, a, lo, hi
  }'
}

bench cavity3d_f32 examples/laminar/cavity3dBenchmark   cavity3d ""              100 1000 3040
bench cavity3d_f16 examples/laminar/cavity3dBenchmark   cavity3d "FEATURES=FP16" 100 1000 6350
bench tgv3d_f32    examples/turbulence/tgv3dBenchmark   tgv3d    ""              300  400 3680
bench tgv3d_f16    examples/turbulence/tgv3dBenchmark   tgv3d    "FEATURES=FP16" 300  400 6900
