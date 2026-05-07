---
name: benchmark-validation
description: Validate solver edits against this repository's convergence behavior and benchmark outputs. Use when changing numerics, boundary conditions, post-processing, grid metrics, MPI reductions, or performance-sensitive code and you need to check `errorU`, `errorT`, `Nu`, `Re`, wall extrema, or the documented `2DRB_ISLBM` regression cases.
---

# Benchmark Validation

Use this skill after meaningful solver edits so validation is deliberate instead of ad hoc.

## Validation ladder

Run the lightest useful validation first, then go deeper only if needed:

1. Syntax check or compile.
2. Short smoke run.
3. Compare against the nearest baseline output.
4. Compare against a documented benchmark if the case matches.

## Repository diagnostics to check

Prioritize the diagnostics already emitted by the core solvers:

- `errorU`
- `errorT`
- `Nu_global`
- `Nu_hot`
- `Nu_cold`
- `Nu_middle`
- `NuVolAvg`
- `ReVolAvg`
- wall-extrema positions
- centerline or mid-plane velocity maxima when relevant

If a change affects only one family of metrics, say that clearly instead of claiming the whole solver is validated.

## Current benchmark anchor

Use `docs/ISLBM_BENCHMARK.md` when validating `2DRB_ISLBM.F90` or closely related 2D side-heated non-uniform-mesh work.

The documented regression anchor there is:

- `Ra = 1e6`
- `Pr = 0.71`
- `Ma = 0.1`
- `41x41` and `81x81` practical laptop-scale checks

Treat those cases as regression checks, especially for thermal metrics. Do not overstate the velocity agreement documented there.

## Comparison rules

- Compare the same geometry and boundary-condition family.
- Compare the same macro configuration.
- Compare the same grid or clearly state the grid change.
- Re-baseline metrics if feature length, `deltaT`, or `heatFluxScale` changed by design.
- For MPI cases, verify the reduction path before trusting any global scalar.

## Common failure modes

- Mix a physics change and a post-processing change, then misread the metric shift.
- Compare side-heated output to a Rayleigh-Benard reference.
- Trust rank-local values from an MPI run as if they were global.
- Treat a successful compile as a validated benchmark result.

## Final reporting rule

When using this skill, end with three explicit statements:

- what was actually checked
- what still remains unverified
- whether a mismatch looks most like a physics, post-processing, or parallel-assembly issue
