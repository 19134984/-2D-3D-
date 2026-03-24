# 2DRB Code Guide

This folder contains the main 2D side-heated cavity codes and a few supporting variants.

## Main source files

- `2DRB.F90`
  Baseline main program in the original project layout.

- `2DRB_ISLBM.F90`
  Non-uniform-mesh / ISLBM-oriented main source used for the current sidewall natural convection work.

- `ChaiPRE2020.F90`
  Chai reference-style source kept here for comparison and verification.

## Subfolders

- `diagnostics/`
  Small diagnostic and sweep variants of the Chai code:
  - `ChaiPRE2020_41.F90`
  - `ChaiPRE2020_81.F90`
  - `ChaiPRE2020_sweep.F90`

- `docs/`
  Notes and validation records:
  - `ISLBM_BENCHMARK.md`

- `build/`
  Local temporary build artifacts for the baseline path. Ignored by git.

- `build_islbm/`
  Local temporary build artifacts for the ISLBM path. Ignored by git.

## Suggested entry points

- If you want the original baseline code, start with `2DRB.F90`.
- If you want the non-uniform mesh / ISLBM version, start with `2DRB_ISLBM.F90`.
- If you want the Chai comparison code, start with `ChaiPRE2020.F90`.
- If you want quick small-grid checks or parameter sweeps, look in `diagnostics/`.
