# Thermal Turbulence / Thermal LBM Literature Download Status

Date: 2026-04-28

Workflow:
- Discovery: web search plus installed `paper-lookup` / `paper-fetch` workflow.
- Download: legal open-access sources only (`paper-fetch` via Semantic Scholar, arXiv, Europe PMC/PMC, plus checked public publisher/repository links).
- Deduplication baseline: existing `pdf/*.pdf`, `pdf/download_status.md`, and `pdf/rb2d_closed_cavity_download_status.md`.
- `UNPAYWALL_EMAIL` was unset, so `paper-fetch` skipped Unpaywall and relied on the other legal OA sources.

## Downloaded successfully

| DOI | Why selected | Local PDF |
|---|---|---|
| `10.1103/RevModPhys.81.503` | Very high-citation RB review; heat transfer and large-scale dynamics baseline | `Ahlers_2008_Heat_transfer_and_large_scale_dynamics_i.pdf` |
| `10.1140/epje/i2012-12058-1` | Good review of turbulent RB convection and open questions | `Chill_2012_New_perspectives_in_turbulent_Rayleigh_B.pdf` |
| `10.1038/s41467-018-04478-0` | Nature Communications; turbulent superstructures in RB | `Pandey_2018_Turbulent_superstructures_in_Rayleigh_B_.pdf` |
| `10.1088/1367-2630/12/7/075022` | Boundary-layer resolution criteria for turbulent thermal convection | `Shishkina_2010_Boundary_layer_structure_in_turbulent_th.pdf` |
| `10.1073/pnas.1922794117` | PNAS; very high Rayleigh number scaling to `Ra = 10^15` | `Iyer_2020_Classical_1_3_scaling_of_convection_hold.pdf` |
| `10.1016/j.jcp.2014.05.041` | Detailed MRT-LBM convection-diffusion model; useful for thermal LBM implementation | `Huang_2014_A_modified_multiple_relaxation_time_latt.pdf` |
| `10.1007/s00366-026-02288-3` | Recent 2026 FVM-LBM coupling for heat transfer simulations | `Zhou_2026_Coupling_finite_volume_lattice_Boltzmann.pdf` |
| `10.1103/PhysRevE.111.035312` | Recent volumetric LBM for thermal particulate flows with conjugate heat transfer | `Zhang_2024_Volumetric_lattice_Boltzmann_method_for_.pdf` |
| `10.1073/pnas.2502972122` | Recent PNAS thermal-plume network paper; detailed thermal turbulence dynamics | `Shevkar_2025_Hierarchical_network_of_thermal_plumes_a.pdf` |

All nine files were checked by magic bytes and begin with `%PDF-`.

## Located but not downloaded automatically

| DOI | Title / reason selected | Status |
|---|---|---|
| `10.1016/j.jcp.2012.11.027` | Boundary conditions for thermal lattice Boltzmann equation method; important thermal LBM boundary paper | No OA PDF found by `paper-fetch` public mode |
| `10.1016/j.jcp.2010.06.037` | Multiple-relaxation-time LBM for convection and anisotropic diffusion; classic detailed MRT-CDE paper | No OA PDF found by `paper-fetch` public mode |
| `10.1016/j.ijthermalsci.2017.04.020` | Cascaded thermal LBM for advection-diffusion and convective heat transfer | No OA PDF found by `paper-fetch` public mode |
| `10.1016/j.compfluid.2018.08.021` | Entropic LBM simulation of thermal convective turbulence | No OA PDF found by `paper-fetch` public mode |
| `10.1016/j.compfluid.2024.106268` | 2024 3D particle-laden thermal convection using LBM | No OA PDF found by `paper-fetch` public mode |
| `10.1017/jfm.2025.243` | JFM 2025 LBM simulations of RB with compressibility-induced non-Oberbeck-Boussinesq effects | Cambridge public PDF link was located, but direct command-line download timed out |
| `10.1016/j.ijheatmasstransfer.2023.124167` | LBM heat transfer in transitional flows with unified single-node curved boundary conditions | No OA PDF found by `paper-fetch` public mode |
| `10.1016/j.tsep.2025.104280` | Recent sidewall oscillation control of RB heat transfer and flow structure | No OA PDF found by `paper-fetch` public mode |
| `10.1016/j.icheatmasstransfer.2024.107525` | Hybrid LBM/immersed-boundary/finite-difference model for thermal fluid-solid interactions | DOI/repository landing pages found, but scripted response was HTML, not PDF |
| `10.1103/PhysRevE.111.055304` | Recent MRT-LBM for conjugate heat transfer | No OA PDF found by `paper-fetch` public mode |

## Notes for follow-up reading

Suggested first pass:
1. `Ahlers_2008_...` and `Chill_2012_...` for broad thermal turbulence/RB context.
2. `Shishkina_2010_...`, `Iyer_2020_...`, and `Shevkar_2025_...` for boundary-layer, high-Ra, and plume-network physics.
3. `Huang_2014_...`, `Zhou_2026_...`, and `Zhang_2024_...` for thermal LBM / heat-transfer implementation ideas.
