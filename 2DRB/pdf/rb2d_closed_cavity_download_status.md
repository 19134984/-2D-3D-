# 2D Closed-Cavity / Square-Cell RB Turbulence Literature Status

Date: 2026-04-24

Workflow used:
- Discovery and metadata checks: paper-lookup workflow with OpenAlex, Crossref, Semantic Scholar, Cambridge Core, PubMed.
- Legal OA download: paper-fetch first; direct verified OA URLs only when paper-fetch exposed a candidate or the publisher/search record showed an open PDF.
- No paywall bypass was used. HTML login pages, bot challenges, and non-PDF responses were rejected.

## Downloaded PDFs

| DOI | Journal | Topic | Local PDF |
|---|---|---|---|
| 10.1063/1.4744988 | Physics of Fluids | 2D RBC flow states versus aspect ratio and Ra | `Poel_2012_Flow_states_in_two_dimensional_Rayleigh_.pdf` |
| 10.1103/PhysRevLett.110.114503 | Physical Review Letters | Flow reversals via vortex reconnections | `Chandra_2012_Flow_reversals_in_turbulent_convection_v.pdf` |
| 10.1103/PhysRevLett.120.144502 | Physical Review Letters | Ultimate-regime transition in 2D RBC | `Zhu_2018_Transition_to_the_Ultimate_Regime_in_Two.pdf` |
| 10.1103/PhysRevE.84.045303 | Physical Review E | Flow structures and heat flux in turbulent RBC | `Poel_2011_Connecting_flow_structures_and_heat_flux.pdf` |
| 10.1063/5.0024408 | Physics of Fluids | Internal flow structure and heat-transfer efficiency | `Xu_2020_Correlation_of_internal_flow_structure_w.pdf` |
| 10.1063/5.0188430 | Physics of Fluids | POD analysis of turbulent RBC flow | `Olesen_2023_Dissipation_based_proper_orthogonal_deco.pdf` |

## High-Priority Papers Located But Not Downloaded Automatically

| DOI | Journal | Why important | Automated download status |
|---|---|---|---|
| 10.1017/jfm.2015.15 | Journal of Fluid Mechanics | 2D square-cell wind reversal at turbulent Ra | OpenAlex marks closed; Cambridge direct PDF request returned HTML |
| 10.1103/PhysRevE.96.023105 | Physical Review E | Velocity/temperature fluctuation statistics in 2D RBC, Ra 1e6-1e10 | OpenAlex marks closed; APS direct request returned 403 |
| 10.1103/PhysRevE.95.013112 | Physical Review E | Precursor model for square-cell wind reversal | OpenAlex marks closed; APS direct request returned 403 |
| 10.1017/jfm.2018.451 | Journal of Fluid Mechanics | 2D tilted-cell flow reversals in rectangular cavities | OpenAlex marks closed; Cambridge direct PDF request returned HTML |
| 10.1017/jfm.2016.181 | Journal of Fluid Mechanics | Quasi-2D geometric confinement experiments | OpenAlex marks closed |
| 10.1017/jfm.2016.647 | Journal of Fluid Mechanics | Reversal cycle in square RB cells | HAL OA candidate found; automated request hit bot/non-PDF page |
| 10.1017/jfm.2019.598 | Journal of Fluid Mechanics | Cessations and reversals in 2D square RB cells | HAL OA candidate found; automated request hit bot/non-PDF page |
| 10.1017/jfm.2020.210 | Journal of Fluid Mechanics | Sidewall temperature control of 2D flow reversal | Cambridge page found; no automated PDF access |
| 10.1017/jfm.2021.58 | Journal of Fluid Mechanics | Stabilizing/destabilizing LSC with sidewall control | Cambridge direct PDF request returned HTML |
| 10.1017/jfm.2022.602 | Journal of Fluid Mechanics | Suppressing reversals by manipulating corner rolls | Cambridge direct PDF request returned HTML |
| 10.1017/jfm.2024.388 | Journal of Fluid Mechanics | Partially isothermal plates, flow reversal and multiple states | Cambridge direct PDF request returned HTML |
| 10.1063/1.5081031 | Chaos | Embedding/POD analysis of reversals in 2D turbulent RBC | HAL OA candidate found; automated request hit bot/non-PDF page |
| 10.1080/14685248.2021.1916023 | Journal of Turbulence | Sidewall control of large-scale flow/reversal | Figshare OA candidate found; automated request hit WAF challenge |

## Notes

- `UNPAYWALL_EMAIL` was not set, so paper-fetch skipped Unpaywall and relied on Semantic Scholar, arXiv, PMC-style routes.
- The best downloaded starting set for code-relevant reading is the PoF 2012 flow-state paper, PRL 2013 vortex-reconnection reversal paper, and PRL 2018 ultimate-regime paper.
- The best not-yet-downloaded square-cell targets are the JFM 2015, JFM 2016, JFM 2019, and JFM 2020-2024 control papers.
