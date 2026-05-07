---
name: case-clone-from-core
description: "Create new solver variants or feature branches by starting from the closest approved core baseline in this repository: `2DRB.F90`, `3DRB.F90`, `3DRBOpenacc.F90`, or `3DRBOpenaccMpi.F90`. Use when the user asks for a new case, port, experimental branch, or derived solver and Codex should choose the right parent file instead of inventing a new layout."
---

# Case Clone From Core

Use this skill when the task is "make a new variant" rather than "edit the existing file in place."

## Parent-selection rule

Choose exactly one primary parent unless the user explicitly asks for a hybrid:

- `2DRB.F90` for 2D cases
- `3DRB.F90` for 3D CPU cases
- `3DRBOpenacc.F90` for 3D GPU cases
- `3DRBOpenaccMpi.F90` for 3D MPI + GPU cases

Follow `AGENTS.md` if there is any doubt.

## Clone workflow

1. State which parent was chosen and why.
2. Copy the smallest viable solver frame from that parent.
3. Preserve naming style, macro layout, output conventions, and post-processing structure unless the new task specifically changes them.
4. List every intentional departure from the parent.
5. Keep reference implementations in `references/` as comparison material, not as the new architecture, unless the user explicitly asks otherwise.

## What to preserve by default

Preserve these repository conventions whenever practical:

- top-level case-selection macros
- main solver phase ordering
- array naming and loop orientation
- logging and settings-file style
- existing `Nu`, `Re`, and convergence output paths

## When to pause and call out risk

Pause before proceeding if the requested change would require:

- mixing two different core parents
- changing the decomposition model
- changing the physical meaning of existing diagnostics
- replacing the established program structure with a brand-new design

Make the departure explicit instead of hiding it inside a clone.

## Common failure modes

- Start from the wrong dimensionality just because a reference file looks closer.
- Pull in features from `references/` that do not match the approved core structure.
- Rename too aggressively and lose the link to the parent solver.
- Change file layout and post-processing at the same time as the physics change.
