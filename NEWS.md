# AIBECS.jl release notes

## v0.15.0 — Package extensions

This release moves a large fraction of AIBECS's heavy dependencies from
`[deps]` to `[weakdeps]` (Julia 1.10+ package extensions), so that
`using AIBECS` no longer transitively pulls plotting / DataFrames /
distribution / NetCDF / JLD2 / MAT / Shapefile / Bijectors stacks.
This significantly cuts time-to-first-plot and time-to-first-load for
users who don't need every feature.

### Breaking changes

You now need to explicitly load a trigger package to enable each
feature group. The relevant `load()` and helper functions emit a
friendly error message naming the required packages until the matching
extension is activated.

| Feature | Required `using ...` |
|---|---|
| Plotting recipes (`plothorizontalslice`, `surfacemap`, `plottransect`, etc.) | `Plots` (or anything that loads `RecipesBase`) |
| Cruise / vertical-profile recipe (`RatioAtStation`) | `Plots` + `Interpolations` |
| Parameter PDF recipes (`PlotParameter`, `PlotParameters`) | `Plots` + `Distributions` |
| `bijector(MyParams)`, `p2λ`, `λ2p`, `mismatch(p)` | `Bijectors` |
| `mismatch(::Distribution, λ)` and prior-based parameter mismatch | `Bijectors` + `Distributions` |
| `table(p)`, `latex(p)` | `DataFrames` |
| `OCIM0.load()`, `OCIM1.load()`, `OCIM2.load()`, `OCCA.load()` | `JLD2` |
| `AeolianSources.load(...)` | `NCDatasets` |
| `ETOPO.load`, `ETOPO.fractiontopo`, `ETOPO.roughnesstopo` | `Distances` + `NCDatasets` |
| `OCIM2_48L.load()` | `MAT` + `NCDatasets` |
| `GroundWaters.load()` | `Shapefile` + `DataFrames` |

### Other behaviour changes

- `Base.show(::AbstractParameters)` now prints a simple per-line
  `name = value (unit)` listing in the default environment. To get the
  prior `DataFrame`-style table, load `using DataFrames` and call
  `table(p)` (or `latex(p)`).
- The `bijector` symbol is no longer exported by AIBECS; use
  `using Bijectors; bijector(MyParams)` (the same function, exported
  by Bijectors.jl). `p2λ`, `λ2p`, and `table` remain exported by
  AIBECS.
- The `@units`, `@prior`, `@flattenable`, `@initial_value`,
  `@description`, `@bounds`, `@logscaled`, `@reference` macros still
  work the same way — `import AIBECS: @units, @prior, ...` and chain
  them with bare names. To use `@prior MyParams ... | LogNormal(0,1) | ...`
  you must now also `using Distributions` explicitly (it is no longer
  pulled transitively by `using AIBECS`).
- Three legacy dead deps were removed: `ImageFiltering`,
  `NearestNeighbors`, `CodecZlib` (JLD2 0.5 pulls
  `ChunkCodecLibZlib` itself).
- `julia` compat bumped to `1.10` (required by the `[extensions]`
  table).
