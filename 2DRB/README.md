# 2DRB Code Guide

This repository now groups solver variants by numerical method and role. Prefer
working inside the matching folder instead of adding new top-level source files.

## Main Folders

- `均匀网格/`
  Uniform-grid 2D/3D solvers. This is the default baseline family for new
  natural-convection, Rayleigh-Benard, OpenMP, OpenACC, and MPI work.

- `ISLBM/`
  Interpolation-supplemented / non-uniform-mesh LBM variants. Use this folder
  for ISLBM development and 2D-to-3D ISLBM alignment work.

- `STLBM/`
  Simplified thermal LBM variants.

- `Xs/`
  Latest LBM-CDE reproduction sources with the `chi` parameter.

- `基准代码/`
  Small reference or benchmark baseline code kept for comparison.

- `后处理/`
  Post-processing programs.

- `运行脚本/`
  Shell and PBS run scripts.

- `references-code/`
  External or paper/reference implementations. Treat these as comparison
  material unless a task explicitly asks to adopt one as the main framework.

- `pdf/`
  The currently retained source papers used by the project.

- `文档/`
  Chinese research reports, derivations, literature evidence, validation
  boundaries, and the merged `最终成果/` deliverables. Start with
  `文档/README.md` for a short guide.

- `tools/`
  Local checking and helper scripts.

## Core Baselines

Use these files as the canonical uniform-grid starting points unless the task
explicitly targets another method family:

- `均匀网格/2DRBOpenmp.F90`
  2D OpenMP baseline.

- `均匀网格/3DRBOpenmp.F90`
  3D OpenMP baseline.

- `均匀网格/3DRBOpenacc.F90`
  3D single-GPU OpenACC baseline.

- `均匀网格/3DRBOpenaccMpi.F90`
  3D MPI + OpenACC baseline.

## Method-Specific Entry Points

- ISLBM 2D CPU: `ISLBM/2DRBOpenmpISLBM.F90`
- ISLBM 3D CPU: `ISLBM/3DRBOpenmpISLBM.F90`
- ISLBM 3D MPI CPU: `ISLBM/3DRBOpenmpMpiISLBM.F90`
- ISLBM 3D GPU: `ISLBM/3DRBOpenaccISLBM.F90`
- ISLBM 3D MPI + GPU: `ISLBM/3DRBOpenaccMpiISLBM.F90`
- STLBM 2D CPU: `STLBM/2DRBOpenmpSTLBM.F90`
- STLBM 2D GPU: `STLBM/2DRBOpenaccSTLBM.F90`
- Flow D2Q9 TRT + thermal D2Q9 BGK LBM-CDE (OpenACC):
  `Xs/D2Q9TRT_D2Q9BGK/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9BGK.F90`
- Flow D2Q9 TRT + fixed-rate thermal D2Q9 TRT LBM-CDE (OpenACC):
  `Xs/D2Q9TRT_D2Q9TRT/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.F90`
- Flow D2Q9 TRT + thermal D2Q5 Luo/Wang TRT (OpenACC):
  `Xs/D2Q9TRT_D2Q5LuoTRT/2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90`

All `Xs` solvers use the `chi_s`-corrected effective flow collision scale:
`(1-chi_s)(tau_even-0.5)(tau_odd-0.5)=3/16`. The former raw/base magic branch
is no longer supported, and no flow-magic compile switch is required.

## Build Helper

The local build helper accepts subfolder paths:

```powershell
powershell -ExecutionPolicy Bypass -File .\build_in_ascii_path.ps1 -SourceFile .\均匀网格\2DRBOpenmp.F90 -SyntaxOnly
powershell -ExecutionPolicy Bypass -File .\build_in_ascii_path.ps1 -SourceFile .\ISLBM\3DRBOpenmpISLBM.F90 -SyntaxOnly
```

Use `-UseMpi` for MPI sources and override `-ParallelFlag` or `-ExtraArgs` only
when the targeted source requires different compiler options.
