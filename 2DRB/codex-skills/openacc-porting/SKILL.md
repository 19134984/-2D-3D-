---
name: openacc-porting
description: Port or review this repository's 3D GPU solvers using the existing OpenACC style from `3DRBOpenacc.F90` and `3DRBOpenaccMpi.F90`. Use when adding `!$acc` directives, creating data regions, handling reductions, synchronizing device and host state, or translating CPU loops to GPU kernels without changing solver physics.
---

# OpenACC Porting

Use this skill when adapting CPU-style kernels to the repository's established OpenACC structure.

## Start from the existing GPU baselines

Treat these files as the source of truth:

- `3DRBOpenacc.F90`
- `3DRBOpenaccMpi.F90`

Use `3DRB.F90` only as the CPU reference for algorithm order and intent.

## Porting workflow

### 1. Preserve algorithm order first

Keep the same high-level sequence as the CPU baseline:

- collision
- streaming
- boundary treatment
- macro reconstruction
- diagnostics and output

Do not rearrange solver phases just to make directives easier to write.

### 2. Define data lifetime before loop directives

Decide which arrays:

- stay resident on device
- need `present(...)`
- need host-to-device initialization
- need device-to-host updates for MPI, logging, or file output

Prefer a coherent data-region strategy over scattered one-off updates.

### 3. Port loops in the repository style

Favor the established pattern:

- `!$acc parallel loop`
- `collapse(...)` for regular nested loops
- explicit `present(...)` clauses for resident arrays
- scalar `reduction(...)` for sums, errors, `Nu`, and `Re`

Keep temporary scalars local to the loop nest and avoid hidden shared state.

### 4. Synchronize host and device at boundaries

Before any host-side operation on device data:

- `MPI` communication
- buffer packing
- restart writes
- diagnostic file output

update the needed slices back to host first, then push changed halos or buffers back to device afterward.

### 5. Validate against the CPU intent

After porting:

- compare outputs against the nearest CPU baseline
- verify `errorU`, `errorT`, `Nu`, and `Re` definitions did not change accidentally
- separate OpenACC issues from physics issues in the final diagnosis

Use `lbm-physics-guard` if the change touches solver meaning, not just execution placement.

## Guardrails

- Do not claim OpenACC verification without an OpenACC-capable compiler.
- Do not rely on a plain Fortran compiler that ignores directives as proof of correctness.
- Do not move host-side file I/O into GPU loops.
- Do not introduce races through reused temporary arrays or missing reductions.
- When MPI is also involved, follow `mpi-domain-decomposition` for halo and reduction rules.
