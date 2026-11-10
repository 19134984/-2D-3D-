#!/usr/bin/env nix
#!nix --extra-experimental-features nix-command --extra-experimental-features flakes develop ../.#gcc-mpi-cuda --command bash

set -euo pipefail

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
export CUDA_VISIBLE_DEVICES=$(( $(nvidia-smi --list-gpus | wc -l) > 1 ? 1 : 0 ))
export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader --id=$CUDA_VISIBLE_DEVICES | grep -o [0-9] | tr -d '\n')

echo "Testing on $(nvidia-smi --query-gpu=name --format=csv,noheader --id=$CUDA_VISIBLE_DEVICES)"

make

ROOT=$PWD

cd examples/fsi/rigidValve2d
make

./rigidValve2d --RESOLUTION 40 --VTK_ENABLED 0

# Correctness: verify the angle(t) transient against the git-tracked reference
# profile (.ci/reference, generated from master) via phase-robust functionals
# (equilibrium offset, amplitude, period) plus an RMS backstop -- see
# compare_valve_profile.py. The overlay+error plot is published as a CI
# artifact. numpy/matplotlib come from an inner nix-shell (not in the build dev
# shell). The last argument is the RMS backstop tolerance (rad).
export MPLCONFIGDIR=$PWD/.mplconfig
nix-shell -p "python3.withPackages(p: [p.numpy p.matplotlib])" --run \
  "python $ROOT/.ci/compare_valve_profile.py \
     $ROOT/.ci/reference/rigidValve2d_n40_angle.csv \
     tmp/gnuplotData/data/valve_n40.csv \
     'gcc-mpi-cuda' 3e-2 valve_profile_n40.png"
