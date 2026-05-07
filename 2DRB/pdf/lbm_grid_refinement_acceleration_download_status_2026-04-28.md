# LBM Grid Refinement / Non-Uniform Grid / Acceleration Download Status

Date: 2026-04-28

Workflow:
- Scope: papers from the previous two search replies on LBM grid refinement, non-uniform grids, AMR/multiresolution grids, algorithmic acceleration, and CPU/GPU/HPC acceleration.
- Discovery/metadata checks: `paper-lookup` style DOI/title verification plus web/source checks where DOI was ambiguous.
- Download: `paper-fetch` public mode first; then `paper-fetch` institutional mode for the public-mode failures. No paywall-bypass sources were used.
- Deduplication baseline: existing `pdf/*.pdf`, `pdf/download_status.md`, `pdf/rb2d_closed_cavity_download_status.md`, and `pdf/thermal_turbulence_lbm_download_status_2026-04-28.md`.
- `UNPAYWALL_EMAIL` was unset, so public mode skipped Unpaywall and relied mainly on Semantic Scholar/arXiv/repository candidates. Institutional mode only let publishers decide access from the current network/IP and used the built-in 1 req/s rate limit.

## DOI manifests

- Full target list: `pdf/lbm_grid_refinement_acceleration_dois_2026-04-28.txt`
- Public-mode failures retried in institutional mode: `pdf/lbm_grid_refinement_acceleration_public_failed_dois_2026-04-28.txt`

## Downloaded successfully

All files below were checked by magic bytes and begin with `%PDF-`.

| DOI | Title / why selected | Local PDF | Source route |
|---|---|---|---|
| `10.1016/j.jcp.2012.03.015` | Advances in multi-domain lattice Boltzmann grid refinement; core multi-domain refinement paper | `Lagrava_2012_Advances_in_multi_domain_lattice_Boltzma.pdf` | Semantic Scholar repository PDF |
| `10.1103/PhysRevE.94.053311` | Grid refinement for entropic lattice Boltzmann models; stability-oriented refinement | `Dorschner_2016_Grid_refinement_for_entropic_lattice_Bol.pdf` | arXiv via Semantic Scholar |
| `10.1016/j.compfluid.2012.02.013` | Performance engineering for the Lattice Boltzmann method on GPGPUs | `Habich_2011_Performance_engineering_for_the_Lattice_.pdf` | arXiv via Semantic Scholar |
| `10.1016/j.cpc.2014.06.003` | Memory transfer optimization for a lattice Boltzmann solver on Kepler architecture NVIDIA GPUs | `Mawson_2013_Memory_transfer_optimization_for_a_latti.pdf` | arXiv via Semantic Scholar |
| `10.1016/j.ijheatmasstransfer.2025.127790` | Interpolation-supplemented LBM simulation of thermal convection on non-uniform meshes | `Xu_2025_Interpolation_supplemented_lattice_Boltz.pdf` | arXiv |
| `10.1016/j.jcp.2020.109645` | Analysis and reduction of spurious noise at grid refinement interfaces | `Astoul_2020_Analysis_and_reduction_of_spurious_noise.pdf` | arXiv via Semantic Scholar |
| `10.1016/j.camwa.2012.05.002` | Comparison of different propagation steps for the lattice Boltzmann method | `Wittmann_2011_Domain_decomposition_and_locality_optimi.pdf` | arXiv via Semantic Scholar |
| `10.1007/s11227-026-08292-0` | TNL-LBM: Scalable LBM implementation based on Template Numerical Library | `Klinkovsk_2026_TNL_LBM_Scalable_lattice_Boltzmann_metho.pdf` | Springer publisher-direct, institutional mode |

## Located but not downloaded automatically

| DOI | Title / reason selected | Automated download status |
|---|---|---|
| `10.1006/jcph.1998.6089` | Grid Refinement for Lattice-BGK Models | No OA PDF found by public/institutional `paper-fetch` routes |
| `10.1103/PhysRevE.67.066707` | Theory and applications of an alternative lattice Boltzmann grid refinement algorithm | No OA PDF found by public/institutional `paper-fetch` routes |
| `10.1103/PhysRevE.89.033310` | Finite-difference LBM with block-structured adaptive-mesh-refinement | No OA PDF found by public/institutional `paper-fetch` routes |
| `10.1016/j.camwa.2015.05.021` | Designing correct fluid hydrodynamics on a rectangular grid using MRT LBM | ScienceDirect publisher-direct candidate found, but scripted request failed / returned non-PDF |
| `10.1109/TPDS.2021.3061895` | Analysis of GPU Data Access Patterns on Complex Geometries for the D3Q19 LBM Algorithm | No OA PDF found by public/institutional `paper-fetch` routes |
| `10.1109/IPDPS57955.2024.00042` | Optimized GPU Implementation of Grid Refinement in LBM | No OA PDF found by public/institutional `paper-fetch` routes |
| `10.3390/fluids6040148` | A Multiple-Grid LBM for Natural Convection under Low and High Prandtl Numbers | MDPI public PDF candidate found, but automated requests returned `403` |
| `10.1007/s10494-025-00689-w` | DNS of turbulent channel/duct flows using interpolation-based LBM | Springer PDF candidate found in institutional mode, but automated request failed |
| `10.1016/j.physa.2005.09.036` | Grid refinement in LBM based on volumetric formulation | ScienceDirect publisher-direct candidate found, but scripted request failed / returned `403` |
| `10.1016/j.compfluid.2013.01.013` | A lattice-Boltzmann method with hierarchically refined meshes | ScienceDirect publisher-direct candidate found, but scripted request failed / returned `403` |
| `10.1016/j.jcp.2013.07.037` | Direct and large-eddy simulation on composite multi-resolution grids by LBM | HAL candidate returned non-PDF; ScienceDirect publisher-direct candidate returned `403` |
| `10.3390/fluids8030103` | Analysis of Hierarchical Grid Refinement Techniques for LBM by Numerical Experiments | MDPI public PDF candidate found, but automated requests returned `403` |
| `10.1006/jcph.1996.0255` | Some Progress in LBM. Part I. Nonuniform Mesh Grids | No OA PDF found by public/institutional `paper-fetch` routes |
| `10.1103/PhysRevE.59.6202` | Finite-volume lattice Boltzmann method | No OA PDF found by public/institutional `paper-fetch` routes |
| `10.4208/cicp.211015.040316a` | FV-LBM for nearly incompressible flows on arbitrary unstructured meshes | No OA PDF found by public/institutional `paper-fetch` routes |
| `10.1016/j.compfluid.2014.11.013` | Numerics of LBM on nonuniform grids: standard LBM and finite-difference LBM | ScienceDirect publisher-direct candidate found, but scripted request failed / returned `403` |
| `10.1016/j.camwa.2026.03.032` | Interpolation-free LBM on non-uniform lattices with second-order recalibration | ScienceDirect publisher-direct candidate found, but scripted request failed / returned `403` |
| `10.1103/PhysRevE.70.066706` | Preconditioned lattice-Boltzmann method for steady flows | No OA PDF found by public/institutional `paper-fetch` routes |
| `10.1016/j.jcp.2009.05.040` | Optimal preconditioning of lattice Boltzmann methods | CSIC / ScienceDirect candidates found, but scripted requests failed |
| `10.1016/j.compfluid.2005.07.020` | Multigrid solution of the steady-state lattice Boltzmann equation | ScienceDirect publisher-direct candidate found, but scripted request failed / returned `403` |
| `10.1103/PhysRevE.101.023309` | Multigrid dual-time-stepping lattice Boltzmann method | HAL candidate returned non-PDF |
| `10.1016/j.camwa.2011.04.012` | Multi-thread implementations of LBM on non-uniform grids for CPUs and GPUs | DOI/SciDirect candidates returned non-PDF or `403` |
| `10.3390/computation5020019` | Esoteric Twist: efficient in-place streaming for LBM | MDPI public PDF candidate found, but automated requests returned `403` |
| `10.3390/computation10060092` | Esoteric Pull and Esoteric Push: in-place streaming schemes for LBM on GPUs | MDPI public PDF candidate found, but automated requests returned `403` |
| `10.1016/j.compfluid.2020.104647` | DNS of turbulent channel flows with mesh-refinement LBM on GPU cluster | ScienceDirect publisher-direct candidate found, but scripted request failed / returned `403` |

## Notes for follow-up

- Highest-priority downloaded papers for immediate reading: `Xu_2025_Interpolation_supplemented_...`, `Lagrava_2012_Advances_...`, `Astoul_2020_Analysis_...`, and `Wittmann_2011_...`.
- Highest-priority still-missing papers: Filippova 1998, Dupuis 2003, Fakhari 2014/2015, Chen 2006, Eitel-Amor 2013, and the MDPI Esoteric Pull/Push + Multiple-Grid papers.
- Several failed items have public landing pages or publisher candidates; browser download or university-library access may work even though command-line retrieval was blocked.
