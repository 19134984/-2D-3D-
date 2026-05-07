# Download status

Date: 2026-04-24

This folder was built from the installed skills workflow:
- `paper-lookup` for discovery and DOI/full-text link checks
- `paper-fetch` for OA PDF dry-run and download attempts

## Downloaded successfully

- [2025_Schupbach_FokkerPlanck_CentralMoment_ThermalConvectiveFlows.pdf](</d:/桌面/代码/代码/热对流/2DRB/pdf/2025_Schupbach_FokkerPlanck_CentralMoment_ThermalConvectiveFlows.pdf>)
  DOI: `10.3390/en18081890`
  Source used: arXiv preprint found via `paper-fetch`

- [2023_Xu_ParticleResolved_ThermalLBM_OpenACC_MultiGPUs.pdf](</d:/桌面/代码/代码/热对流/2DRB/pdf/2023_Xu_ParticleResolved_ThermalLBM_OpenACC_MultiGPUs.pdf>)
  DOI: `10.1016/j.ijheatmasstransfer.2023.124758`
  Source used: arXiv preprint found via `paper-fetch`

- [2018_Molla_GPU_Accelerated_MRT_LBM_ConvectiveFlows_PorousMedia.pdf](</d:/桌面/代码/代码/热对流/2DRB/pdf/2018_Molla_GPU_Accelerated_MRT_LBM_ConvectiveFlows_PorousMedia.pdf>)
  DOI: `10.3389/fmech.2018.00015`
  Source used: Frontiers OA PDF found via `paper-fetch`

- [2017_Fei_CentralMoments_Based_LBM_ThermalFlows.pdf](</d:/桌面/代码/代码/热对流/2DRB/pdf/2017_Fei_CentralMoments_Based_LBM_ThermalFlows.pdf>)
  DOI: `10.1016/j.ijheatmasstransfer.2017.12.052`
  Source used: arXiv preprint found via `paper-fetch`

- [Allen_2016_Moment_based_boundary_conditions_for_lat.pdf](</d:/桌面/代码/代码/热对流/2DRB/pdf/Allen_2016_Moment_based_boundary_conditions_for_lat.pdf>)
  DOI: `10.1504/PCFD.2016.077296`
  Source used: repository PDF downloaded successfully by `paper-fetch`

## Not downloaded in this run

- `10.1016/j.icheatmasstransfer.2024.107653`
  Title: A new Neumann boundary condition scheme for the thermal lattice Boltzmann method
  Status: `paper-fetch` public mode did not find an OA PDF

- `10.3390/en19041005`
  Title: An Improved Lattice Boltzmann Method for Simulating High-Conductivity-Ratio Conjugate Heat Transfer
  Status: legal direct PDF URL was identified via `paper-lookup`, but automated download hit `403 Forbidden`

- `10.1002/cpe.70518`
  Title: A Python/Fortran Implementation of the Lattice-Boltzmann Kernel on Multiple GPU Using the OpenACC Framework
  Status: Crossref reported a PDF link, but the publisher blocked scripted access with `403`

- `10.1016/j.procs.2025.08.232`
  Title: Multi-GPU porting of a phase-change cascaded lattice Boltzmann method for three-dimensional pool boiling simulations
  Status: no OA PDF found in current public-sources run

- `10.1109/access.2020.2971546`
  Title: Lattice Boltzmann Method for Fluid-Thermal Systems: Status, Hotspots, Trends and Outlook
  Status: Crossref reported an IEEE PDF link, but scripted access failed

- `10.3390/en16217229`
  Title: Multiple-Relaxation-Time Lattice Boltzmann Simulation of Soret and Dufour Effects on the Thermosolutal Natural Convection of a Nanofluid in a U-Shaped Porous Enclosure
  Status: legal direct PDF URL was identified via `paper-lookup`, but automated download hit `403 Forbidden`

## Notes

- `asta-skill` was checked first, but this session does not currently have `ASTA_API_KEY` or a registered Asta MCP server, so the live search path used `paper-lookup` plus `paper-fetch`.
- `Sci-Hub` or other paywall-bypass sources were not used.
