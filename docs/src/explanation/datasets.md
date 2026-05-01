# Datasets

AIBECS.jl ships a handful of small pedagogical circulations *inside* the
package (built from scratch in pure Julia, no download needed) and downloads
larger data products on demand: ocean circulation matrices, dust deposition
fields, topography, and a few others. This page catalogs everything the
package can build or fetch.

## Why downloads are deferred to first use

The downloaded datasets are large (the OCIM2 transport matrices are ~29 MB
each, OCIM2_48L is ~554 MB, ETOPO is ~400 MB compressed), and bundling them
would bloat every install for users who only need a subset. AIBECS instead
uses [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl) to declare each
dataset as a *data dependency*: a named record with a URL, a checksum, and
a citation.

The first time you call e.g. `OCIM2.load()`, DataDeps:

1. Checks the local DataDeps cache (`~/.julia/datadeps/AIBECS-OCIM2_CTL_He/` by default) for a copy.
2. If absent, prompts you to accept the licence and citation (auto-accepted on CI when `ENV["DATADEPS_ALWAYS_ACCEPT"] = true`), then downloads the file from the source listed below.
3. Verifies the recorded checksum, refusing to use a corrupted file.
4. Caches the file so subsequent loads are instant and offline-friendly.

This pattern (described in [White et al., 2019](https://doi.org/10.5334/jors.244))
is the same mechanism used by `MLDatasets.jl`, `WordNet.jl`, and other
data-heavy Julia packages. The benefit for AIBECS users: the package stays
small and pure-code, the data lives at a stable URL with a citation, and any
change to the upstream file is caught by the checksum mismatch instead of
silently changing model output.

The toy / pedagogical circulations are different: they are built in memory
from a few constants and the helpers in `CirculationGeneration.jl`, so they
add no install-time cost and need no network access.

## Ocean circulations

### Built from scratch (bundled in AIBECS)

These pedagogical circulations are constructed each time `Module.load()` is
called, using `OceanGrids` plus the `T_advection` / `T_diffusion` helpers in
[src/CirculationGeneration.jl](https://github.com/JuliaOcean/AIBECS.jl/blob/main/src/CirculationGeneration.jl). Nothing is
downloaded or cached.

| Module | Layout | Grid size | Citation | Source | Size (MB) |
|---|---|---|---|---|---|
| `TwoBoxModel` | surface + deep | 1×1×2 | [Sarmiento & Gruber (2006)](https://press.princeton.edu/books/hardcover/9780691017075/ocean-biogeochemical-dynamics) | bundled | 0 |
| `Archer_etal_2000` | 3-box (HL surface, LL surface, deep) on a 6-cell grid | 2×1×3 | [Archer et al. (2000)](https://doi.org/10.1029/1999GB001216) | bundled | 0 |
| `Primeau_2x2x2` | shoebox (5 wet boxes, 3 dry) | 2×2×2 | [Primeau, Intro2TransportOperators](https://github.com/fprimeau/BIOGEOCHEM_TEACHING/blob/master/Intro2TransportOperators.ipynb) | bundled | 0 |
| `Haine_and_Hall_2025` | 9-box (3 latitudes × 3 depths) | 3×1×3 | [Haine & Hall (2002)](https://doi.org/10.1175/1520-0485(2002)032%3C1932:AGTTWM%3E2.0.CO;2); [Haine et al. (2025)](https://doi.org/10.1029/2024MS004637) | bundled | 0 |

### Downloaded (data-based circulations)

Five families of transport matrices and grids are available. All files are
JLD2 (or `.tar.gz` for OCIM2_48L) generated from the upstream MATLAB
distributions by [briochemc/OceanCirculations](https://github.com/briochemc/OceanCirculations).

| Module | Variant | Grid size | Citation | Source | Size (MB) |
|---|---|---|---|---|---|
| `OCIM0` | (single) | 180×90×24 | [DeVries & Primeau (2011)](https://doi.org/10.1175/JPO-D-10-05011.1); [Primeau et al. (2013)](https://doi.org/10.1002/jgrc.20181) | [link](https://figshare.com/ndownloader/files/28336086) | 11 |
| `OCIM1` | CTL | 180×91×24 | [DeVries (2014)](https://doi.org/10.1002/2013GB004739) | [link](https://ndownloader.figshare.com/files/28335882) | 28 |
| `OCIM2` | CTL_He (default) | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336284) | 29 |
| `OCIM2` | CTL_noHe | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336299) | 29 |
| `OCIM2` | KiHIGH_He | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336302) | 29 |
| `OCIM2` | KiHIGH_noHe | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336311) | 28 |
| `OCIM2` | KiLOW_He | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336317) | 29 |
| `OCIM2` | KiLOW_noHe | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336323) | 29 |
| `OCIM2` | KvHIGH_He | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336326) | 29 |
| `OCIM2` | KvHIGH_KiHIGH_noHe | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336329) | 28 |
| `OCIM2` | KvHIGH_KiLOW_He | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336332) | 29 |
| `OCIM2` | KvHIGH_KiLOW_noHe | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336341) | 29 |
| `OCIM2` | KvHIGH_noHe | 180×91×24 | [DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716) | [link](https://ndownloader.figshare.com/files/28336353) | 29 |
| `OCIM2_48L` | base | 180×91×48 | [Holzer, DeVries & de Lavergne (2021)](https://doi.org/10.1038/s41467-021-24648-x) | [link](https://ndownloader.figshare.com/files/28468077) | 554 |
| `OCCA` | (single) | 180×80×10 | [Forget (2010)](https://doi.org/10.1175/2009JPO4043.1) | [link](https://ndownloader.figshare.com/files/28336173) | 5 |

Loading any OCIM matrix requires `using JLD2` (which activates the
`AIBECSJLD2Ext` extension); loading `OCIM2_48L` additionally requires
`using MAT, NCDatasets`.

## Other datasets

Aeolian deposition, topography, river discharge, groundwater discharge, and
the AWESOME-OCIM toolbox. These are hosted by their original maintainers
(except Chien, which lives on Zenodo).

| Module | Variant | What it is | Citation | Source | Size (MB) |
|---|---|---|---|---|---|
| `AeolianSources` | Chien (default) | 2°×2° seasonal aerosol deposition (fires, biofuels, dust, sea salt, biogenics, volcanoes, fossil fuels) | [Chien et al. (2016)](https://doi.org/10.1002/2015GB005334) | [link](https://zenodo.org/records/17869479/files/post.aerosols.2x2.seasonal.nc) | 1 |
| `AeolianSources` | Kok | annual dust deposition by source region (DustCOMM) | [Kok et al. (2021)](https://doi.org/10.5194/acp-21-8169-2021) | [link](https://research.aos.ucla.edu/dustcomm/K21b/DustCOMM_source_region_wetdep_annual_PM20_abs.nc) | 6 |
| `ETOPO` | bedrock | 1-arc-min global relief, bedrock surface (compressed) | [Amante & Eakins (2009)](https://doi.org/10.7289/V5C8276M) | [link](https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gdal.grd.gz) | 402 |
| `ETOPO` | ice | 1-arc-min global relief, ice-surface (compressed) | [Amante & Eakins (2009)](https://doi.org/10.7289/V5C8276M) | [link](https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gdal.grd.gz) | 395 |
| `GroundWaters` | (single) | coastal fresh-groundwater discharge shapefile | [Luijendijk et al. (2020)](https://doi.org/10.1038/s41467-020-15064-8); [PANGAEA dataset](https://doi.org/10.1594/PANGAEA.907641) | [link](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-15064-8/MediaObjects/41467_2020_15064_MOESM8_ESM.zip) | 32 |
| `AO` | master | AWESOME-OCIM MATLAB toolbox (source archive) | [John et al. (2020)](https://doi.org/10.1016/j.chemgeo.2019.119403) | [link](https://github.com/hengdiliang/AWESOME-OCIM/archive/refs/heads/master.zip) | 123 |

Loading these typically requires the matching extension dependencies:
`AeolianSources` needs `NCDatasets`, `ETOPO` needs `Distances, NCDatasets`,
`GroundWaters` needs `Shapefile, DataFrames`.

## Migration to Zenodo (in progress)

The OCIM matrices are currently hosted on FigShare (the table above shows
the live URLs). FigShare's CDN occasionally fails during CI, so these files
are being migrated to Zenodo for reliability. Once the migration is
complete this page will be updated with the new Zenodo URLs and DOIs; the
FigShare links above will remain functional but will no longer be the
canonical source.

## Cache layout

DataDeps stores files under `~/.julia/datadeps/<DataDepName>/` by default.
You can override the location by setting `ENV["DATADEPS_LOAD_PATH"]` before
loading AIBECS. The DataDep names this package registers are:

- `AIBECS-OCIM0.1`, `AIBECS-OCIM1_CTL`, `AIBECS-OCIM2_<variant>` (one per OCIM2 variant), `AIBECS-OCIM2_48L`, `AIBECS-OCCA`
- `AIBECS-Chien_etal_2016`, `AIBECS-Kok_etal_2021`
- `ETOPO_bedrock`, `ETOPO_ice`
- `groundwater_discharge`
- `AWESOME-OCIM`

To force a re-download (e.g. after a checksum mismatch or to test a URL
change), delete the corresponding directory and call `Module.load()` again.
