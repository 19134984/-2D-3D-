---
name: fortran-build-run
description: Compile, syntax-check, and smoke-test local Fortran sources in this repository, especially `2DRB.F90`, `3DRB.F90`, `3DRBOpenacc.F90`, `3DRBOpenaccMpi.F90`, and close variants. Use when Codex needs to choose the right baseline, pick `gfortran` vs `mpifort`, run `build_in_ascii_path.ps1`, do a safe syntax-only check first, or explain what could not be verified because the required compiler or runtime is unavailable.
---

# Fortran Build Run

Use this skill to verify Fortran edits in this repository without inventing a new build flow.

## Baseline-first rule

Choose the closest approved core file before compiling or running:

- `2DRB.F90` for 2D work
- `3DRB.F90` for 3D CPU work
- `3DRBOpenacc.F90` for 3D OpenACC work
- `3DRBOpenaccMpi.F90` for 3D MPI + OpenACC work

Treat reference codes as comparison material unless the user explicitly wants to build them.

## Default workflow

1. Inspect the source for `use mpi`, `!$acc`, and case-selection macros.
2. Prefer a syntax-only pass before a full build.
3. Prefer the repository build helper on Windows so the source is copied to an ASCII-safe temp path.
4. Run a smoke test only when the binary and case size make it practical.
5. Report exactly what was verified and what was not.

## Repository build commands

Prefer `build_in_ascii_path.ps1` for CPU or MPI-capable builds on Windows:

- Syntax check CPU baseline:
  - `powershell -ExecutionPolicy Bypass -File .\build_in_ascii_path.ps1 -SourceFile 3DRB.F90 -SyntaxOnly`
- Compile CPU baseline:
  - `powershell -ExecutionPolicy Bypass -File .\build_in_ascii_path.ps1 -SourceFile 3DRB.F90`
- Compile an MPI source that the local MPI wrapper can actually build:
  - `powershell -ExecutionPolicy Bypass -File .\build_in_ascii_path.ps1 -SourceFile <file> -UseMpi`

Use `run.sh` only when a matching executable already exists and the environment is Linux-like.

## Compiler choice rules

- Use `gfortran` through the build helper for plain CPU baselines.
- Use `mpifort` through the build helper only when the target source and wrapper support the needed modules.
- Do not claim OpenACC compile verification unless an OpenACC-capable compiler such as `nvfortran` or `nvc` is actually available.
- Do not treat OpenACC directives as "verified" just because a non-OpenACC compiler ignored them as comments.

## Verification depth

Choose the lightest verification that still answers the task:

- Syntax-only for quick regression checks
- Full compile for code-generation confidence
- Short runtime smoke test for obvious crashes or startup failures
- Longer benchmark runs only when the user asked for validation or performance data

## Guardrails

- Rebuild after changing preprocessor macros, not just Fortran statements.
- Do not trust stale executables if the source, macros, or compiler flags changed.
- Say explicitly when OpenACC or MPI runtime verification was not possible.
- Keep build and run notes short but concrete so later work can reuse them.
