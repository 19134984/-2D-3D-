---
name: lbm-physics-guard
description: Review or implement LBM solver changes while preserving physical and numerical consistency in this repository. Use when editing collision, streaming, equilibrium, source terms, relaxation times, `EnableUseG`, `EnableLegacyThermalScheme`, boundary conditions, feature-length logic, or Nu/Re post-processing in the 2D and 3D core solvers.
---

# LBM Physics Guard

Use this skill when a code change touches the mathematics or physical interpretation of the solver, not just formatting or comments.

## Start from the approved baselines

Check the nearest core solver first:

- `2DRB.F90`
- `3DRB.F90`
- `3DRBOpenacc.F90`
- `3DRBOpenaccMpi.F90`

Treat these files as the source of truth for solver structure, case macros, and post-processing layout.

## Review checklist

### 1. Keep the case definition consistent

When switching or editing `SideHeatedCell` and `RayleighBenardCell`, update the coupled pieces together:

- wall-temperature and velocity boundary macros
- periodic vs adiabatic wall choices
- buoyancy direction and sign convention
- feature-length logic
- `deltaT`, `heatFluxScale`, and any derived comparison scales
- the matching `Nu` and `Re` post-processing routines

If feature length or heat-flux scaling changes, read `docs/CASE_FEATURE_LENGTH_UPDATE_CN.md`.

### 2. Keep lattice-model pieces aligned

When editing D2Q or D3Q kernels, keep the following families consistent:

- lattice weights
- equilibrium definitions
- MRT or relaxation parameters
- forcing or buoyancy source terms
- macro reconstruction formulas

Do not update only one piece of a coupled model.

### 3. Keep dimensionless parameters coupled

When touching `Rayleigh`, `Prandtl`, `Mach`, viscosity, diffusivity, `tauf`, or `gBeta1`, verify that all dependent formulas still agree with the intended physics.

### 4. Keep boundary treatment and diagnostics coupled

When editing bounce-back, isothermal, adiabatic, or periodic logic:

- verify the distribution update itself
- verify the macro field on the affected cells
- verify the downstream `Nu`, `Re`, or extrema calculations that depend on those cells

### 5. Keep convergence metrics meaningful

Preserve the interpretation of:

- `errorU`
- `errorT`
- `Nu_global`
- `Nu_hot`
- `Nu_cold`
- `Nu_middle`
- `NuVolAvg`
- `ReVolAvg`

If a metric changes by design, say so explicitly instead of presenting it as a regression.

## Common failure modes

- Flip a case macro without updating the matching post-processing path.
- Change a heat-flux definition without updating `heatFluxScale`.
- Port a 2D formula into 3D without checking the wall-normal direction.
- Mix `EnableUseG` and `EnableLegacyThermalScheme` changes without checking coupled temperature logic.
- Modify a source term and forget the reference temperature or sign convention.

## Output expectation

When using this skill, explain:

- what physical invariant or numerical assumption changed
- which coupled sections were checked
- whether the result still matches the nearest approved baseline in intent

Use `benchmark-validation` after significant physics-facing edits.
