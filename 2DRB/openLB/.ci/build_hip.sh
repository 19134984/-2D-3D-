#!/usr/bin/env nix
#!nix --extra-experimental-features nix-command --extra-experimental-features flakes develop ../.#hip --command bash

# Compile-only coverage for the GPU_HIP backend targeting AMD gfx via ROCm/hipcc.
# nixpkgs ROCm (rocmPackages.clr) is AMD-only, so this builds device code for the
# gfx arch (default gfx90a, override with OLB_HIP_ARCH) without requiring an AMD
# GPU on the runner. A curated subset of examples exercises the major dynamics,
# coupling and field-storage paths.

set -euo pipefail

make -j4

make -j4 \
     examples/laminar/cavity3d/ \
     examples/laminar/cavity3dBenchmark/ \
     examples/laminar/poiseuille3d/ \
     examples/laminar/bstep3d/ \
     examples/laminar/cylinder3d/ \
     examples/turbulence/nozzle3d/ \
     examples/turbulence/tgv3d/ \
     examples/thermal/rayleighBenard3d/ \
     examples/advectionDiffusionReaction/advectionDiffusion3d/ \
     examples/multiComponent/phaseSeparation3d/
