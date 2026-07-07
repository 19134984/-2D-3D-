# AGENTS.md

## Core Code Baselines

The canonical core source files for this repository are:

- `均匀网格/2DRBOpenmp.F90`
- `均匀网格/3DRBOpenmp.F90`
- `均匀网格/3DRBOpenacc.F90`
- `均匀网格/3DRBOpenaccMpi.F90`

Unless the user explicitly says otherwise, any newly generated code, variant,
port, optimization, refactor, diagnostic case, or bug fix should use these
four files as the main framework.

## Baseline Selection Rule

Choose the closest baseline before generating or editing code:

- `均匀网格/2DRBOpenmp.F90` for 2D uniform-grid natural-convection / Rayleigh-Benard style work
- `均匀网格/3DRBOpenmp.F90` for 3D uniform-grid CPU baseline work
- `均匀网格/3DRBOpenacc.F90` for 3D uniform-grid OpenACC / GPU work
- `均匀网格/3DRBOpenaccMpi.F90` for 3D uniform-grid MPI + OpenACC work

For method-specific variants, start from the closest organized folder:

- `ISLBM/` for interpolation-supplemented / non-uniform-mesh LBM variants
- `STLBM/` for simplified thermal LBM variants
- `Xs/` for the latest LBM-CDE reproduction with the `chi` parameter

## Default Implementation Guidance

- Prefer incremental changes on top of the nearest baseline instead of creating
  a brand-new architecture.
- Preserve the overall program structure, naming style, macro switches,
  physical-model assumptions, and boundary-condition organization from the
  selected baseline whenever practical.
- Treat `references-code/` and other auxiliary codes as comparison material unless
  the user explicitly requests adopting them as the new main framework.
- If a requested change would require a large structural departure from one of
  the four core baselines, call out that departure before proceeding.

## Language

- Use Simplified Chinese for user-facing explanations, progress updates, and
  final responses by default.
- Keep code, commands, file paths, variable names, compiler messages, and API
  identifiers in their original language unless the user explicitly asks for a
  translation.

## Priority

If there is any conflict between older notes and current work, follow this
priority order unless the user provides a different instruction:

1. Explicit user request in the current conversation
2. This `AGENTS.md`
3. Other repository notes or reference implementations
