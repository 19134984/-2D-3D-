#!/usr/bin/env nix
#!nix --extra-experimental-features nix-command --extra-experimental-features flakes develop ../.#gcc-mpi-cuda-vtk --command bash

set -euo pipefail

export CUDA_ARCH=86

cd external/voreen/backend
make
