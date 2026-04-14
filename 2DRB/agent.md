# agent.md

This file mirrors the active repository instructions in `AGENTS.md`.

## Core Code Baselines

The canonical core source files for this repository are:

- `2DRB.F90`
- `3DRB.F90`
- `3DRBOpenacc.F90`
- `3DRBOpenaccMpi.F90`

Unless the user explicitly says otherwise, any newly generated code, variant,
port, optimization, refactor, diagnostic case, or bug fix should use these
four files as the main framework.

## Baseline Selection Rule

Choose the closest baseline before generating or editing code:

- `2DRB.F90` for 2D natural-convection / Rayleigh-Benard style work
- `3DRB.F90` for 3D CPU baseline work
- `3DRBOpenacc.F90` for 3D OpenACC / GPU work
- `3DRBOpenaccMpi.F90` for 3D MPI + OpenACC work

## Default Implementation Guidance

- Prefer incremental changes on top of the nearest baseline instead of creating
  a brand-new architecture.
- Preserve the overall program structure, naming style, macro switches,
  physical-model assumptions, and boundary-condition organization from the
  selected baseline whenever practical.
- Treat `references/` and other auxiliary codes as comparison material unless
  the user explicitly requests adopting them as the new main framework.
- If a requested change would require a large structural departure from one of
  the four core baselines, call out that departure before proceeding.

## Priority

If there is any conflict between older notes and current work, follow this
priority order unless the user provides a different instruction:

1. Explicit user request in the current conversation
2. `AGENTS.md`
3. Other repository notes or reference implementations
