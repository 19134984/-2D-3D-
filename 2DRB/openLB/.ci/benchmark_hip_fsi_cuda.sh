#!/usr/bin/env nix
#!nix --extra-experimental-features nix-command --extra-experimental-features flakes develop ../.#hip-cuda --command bash

# FSI runtime coverage for the GPU_HIP backend targeting NVIDIA
# (HIP_PLATFORM=nvidia). Mirrors benchmark_fsi_cuda.sh but routes the lattice
# through HIP, so the HIP-on-CUDA path gets the same FSI runtime check as native
# CUDA. Sourced entirely from nixpkgs via the hip-cuda dev shell (which also puts
# the real CUDA driver on LD_LIBRARY_PATH), so it runs on the NVIDIA has-gpu
# runner without a system ROCm install and without per-job driver exports.

set -euo pipefail

export CUDA_VISIBLE_DEVICES=$(( $(nvidia-smi --list-gpus | wc -l) > 1 ? 1 : 0 ))

echo "Testing on $(nvidia-smi --query-gpu=name --format=csv,noheader --id=$CUDA_VISIBLE_DEVICES) (GPU_HIP -> nvidia, sm_$CUDA_ARCH)"

make -j4

ROOT=$PWD

pushd examples/fsi/rigidValve2d
make
./rigidValve2d --RESOLUTION 40 --VTK_ENABLED 0

# Correctness: verify the angle(t) transient against the same git-tracked
# reference profile (.ci/reference, from master) as the native FSI job. This
# single-mode HIP build integrates the rigid body with nvcc-compiled host code,
# whose floating-point contraction differs from the mixed CUDA job's gcc host
# pass and phase-shifts the early swing by ~0.04 s (a benign, re-synchronizing
# offset). We therefore gate on phase-robust functionals (equilibrium offset,
# amplitude, period) plus an RMS backstop rather than pointwise max|error| --
# see compare_valve_profile.py. The last argument is the RMS backstop tol (rad).
export MPLCONFIGDIR=$PWD/.mplconfig
nix-shell -p "python3.withPackages(p: [p.numpy p.matplotlib])" --run \
  "python $ROOT/.ci/compare_valve_profile.py \
     $ROOT/.ci/reference/rigidValve2d_n40_angle.csv \
     tmp/gnuplotData/data/valve_n40.csv \
     'hip-cuda' 3e-2 valve_profile_n40.png"
popd
