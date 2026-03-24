# ISLBM Benchmark Notes

This note records the current small-grid validation status of [2DRB_ISLBM.F90](D:\桌面\codex\-2D-3D--main\2DRB\2DRB_ISLBM.F90).

## Case Setup

- Problem: 2D side-heated square cavity
- Reference paper: [Interpolation-supplemented lattice Boltzmann simulation of thermal convection on non-uniform meshes](https://xuaoxiqi.github.io/publications/043_IJHMT2026_Interpolation.pdf)
- Parameters: `Ra = 1e6`, `Pr = 0.71`, `Ma = 0.1`
- Mesh mode: `erf` non-uniform mesh
- Stretch parameter: `a = 1.5`
- Steady convergence thresholds: `epsU = 1e-9`, `epsT = 1e-9`

## Reference Values

Paper Table 2 gives the following reference values for the `Ra = 1e6` side-heated cavity:

| Source | `<Nu>` | `Nu_hot` | `Nu_middle` | `<Re>` |
| --- | ---: | ---: | ---: | ---: |
| Paper reference | 8.8186 | 8.8206 | 8.8186 | 99.3542 |

## Laptop-Scale Results

| Grid | `L0` | `Nu_global` | `Nu_hot` | `Nu_middle` | `Re` | CPU time (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `41x41` | 59.3241 | 8.56906279 | 9.38065001 | 8.48139828 | 72.47660064 | 49.5469 |
| `81x81` | 118.2970 | 8.83872670 | 9.02850783 | 8.83154644 | 73.38021934 | 347.2969 |

## Relative Error vs Paper

| Grid | `<Nu>` err. | `Nu_hot` err. | `Nu_middle` err. | `<Re>` err. |
| --- | ---: | ---: | ---: | ---: |
| `41x41` | `-2.83%` | `+6.35%` | `-3.82%` | `-27.05%` |
| `81x81` | `+0.23%` | `+2.36%` | `+0.15%` | `-26.14%` |

## Interpretation

- `81x81` is already good enough for laptop-scale regression on the thermal side: `Nu_global` and `Nu_middle` are both within about `0.25%` of the paper reference.
- `Nu_hot` is still slightly high on coarse grids, but the error is already down to about `2.4%` on `81x81`.
- `<Re>` remains clearly under-predicted on `41x41` and `81x81`. On this laptop, these two grids should be treated as smoke/regression benchmarks, not as full paper-grade velocity validation.
- For this machine, `41x41` and `81x81` are practical default checks. `257x257` is possible in principle, but is too heavy for routine use on a laptop.

## Saved Run Logs

The numbers above were taken from these local run directories:

- `C:\codextmp\islbm41bench_build_v3`
- `C:\codextmp\islbm81bench_build_v3`
