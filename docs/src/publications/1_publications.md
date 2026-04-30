# [Publications and software using AIBECS](@id publications)

This page collects papers and software that use AIBECS, plus related
projects in the same ecosystem. It is community-driven — if your work
uses AIBECS and is missing here, please open a pull request (see
[How to add your work](#how-to-add-your-work) below).

## Software citation

This code is © Benoît Pasquier (2026) and contributors, and it is made available under the MIT license enclosed with the software.

Over and above the legal restrictions imposed by this license, if you use AIBECS for an academic publication then you are obliged to provide proper attribution. This can be to this code directly,

> Benoît Pasquier, François W. Primeau, and Seth G. John (2026). **AIBECS.jl: A tool for exploring global marine biogeochemical cycles**. _Zenodo_. doi: [10.5281/zenodo.2864051](https://doi.org/10.5281/zenodo.2864051).

or to the paper that describes it,

> Pasquier, B., Primeau, F. W., and John, S. G. (2022). **AIBECS.jl: A tool for exploring global marine biogeochemical cycles**. _Journal of Open Source Software_, 7(69), 3814. doi: [10.21105/joss.03814](https://doi.org/10.21105/joss.03814).

but, ideally, both. BibTeX entries live in [`CITATION.bib`](https://github.com/JuliaOcean/AIBECS.jl/blob/main/CITATION.bib) under the keys `AIBECS.jl` and `Pasquier_etal_JOSS_2022`. If you also use data provided by AIBECS (such as the OCIM circulations), please cite those datasets as well.

## Papers using AIBECS

This list is curated by the maintainer; please add your own work via PR.

- **Du, J., Haley, B. A., McManus, J., Blaser, P., Rickli, J., & Vance,
  D. (2025).** Abyssal seafloor as a key driver of ocean trace-metal
  biogeochemical cycles. *Nature*, 642(8068), 620–627.
  [doi:10.1038/s41586-025-09038-3](https://doi.org/10.1038/s41586-025-09038-3)
- **Liang, Z., Letscher, R. T., & Knapp, A. N. (2025).** Oligotrophic
  ocean new production supported by lateral transport of dissolved
  organic nutrients. *Global Biogeochemical Cycles*, 39(6).
  [doi:10.1029/2024GB008345](https://doi.org/10.1029/2024GB008345)
- **John, S. G., Liang, H., Pasquier, B., Holzer, M., & Silva, S.
  (2024).** Biogeochemical fluxes of nickel in the global oceans
  inferred from a diagnostic model. *Global Biogeochemical Cycles*,
  38(5).
  [doi:10.1029/2023GB008018](https://doi.org/10.1029/2023GB008018)
- **Du, J. (2023).** SedTrace 1.0: a Julia-based framework for
  generating and running reactive-transport models of marine sediment
  diagenesis specializing in trace elements and isotopes.
  *Geoscientific Model Development*, 16(20), 5865–5894.
  [doi:10.5194/gmd-16-5865-2023](https://doi.org/10.5194/gmd-16-5865-2023)
- **Pasquier, B., Hines, S. K. V., Liang, H., Wu, Y., Goldstein, S. L.,
  & John, S. G. (2022).** GNOM v1.0: an optimized steady-state model
  of the modern marine neodymium cycle. *Geoscientific Model
  Development*, 15(11), 4625–4656.
  [doi:10.5194/gmd-15-4625-2022](https://doi.org/10.5194/gmd-15-4625-2022)
- **Pasquier, B., Primeau, F. W., & John, S. G. (2022).** AIBECS.jl:
  A tool for exploring global marine biogeochemical cycles.
  *Journal of Open Source Software*, 7(69), 3814.
  [doi:10.21105/joss.03814](https://doi.org/10.21105/joss.03814)

## Related software

- [GNOM](https://github.com/MTEL-USC/GNOM) — optimised steady-state
  model of the modern marine neodymium cycle, built on AIBECS. Software
  DOI: [10.5281/zenodo.5651295](https://doi.org/10.5281/zenodo.5651295).
- [OceanCirculations.jl](https://github.com/briochemc/OceanCirculations) —
  collection of ocean transport matrices distributed via DataDeps; in
  the near future, this and `OceanGrids.jl` may be bundled directly into
  AIBECS rather than living in separate repos.
- [OceanGrids.jl](https://github.com/briochemc/OceanGrids) — grid
  abstractions used throughout AIBECS; same near-future-bundling note.
- [OceanGreensFunctionMethods.jl](https://github.com/ggebbie/OceanGreensFunctionMethods.jl) —
  Julia implementation of Green's-function tracer-pathway methods,
  including the pedagogical 9-box model surfaced in AIBECS as
  `Haine_and_Hall_2025`.

## How to add your work

If your published paper or software uses AIBECS:

1. Open a pull request adding an entry to the relevant list above
   (alphabetised within section by year, newest first).
2. Append the corresponding BibTeX entry to
   [`CITATION.bib`](https://github.com/JuliaOcean/AIBECS.jl/blob/main/CITATION.bib).
   Use a descriptive citekey (e.g. `Author_etal_Journal_Year`) and
   place the entry under the *Papers using AIBECS* section.
3. Briefly describe the role AIBECS played in your work in the PR
   description so the maintainer can sanity-check the categorisation.
