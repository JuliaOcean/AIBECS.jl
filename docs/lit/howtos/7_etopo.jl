#---------------------------------------------------------
# # [Subgrid topography from ETOPO](@id etopo)
#---------------------------------------------------------

# This how-to shows how to obtain subgrid-scale seafloor information
# from the [ETOPO1](https://www.ngdc.noaa.gov/mgg/global/global.html)
# 1-arc-minute global relief model and bin it onto an AIBECS grid.
# Two diagnostics are exposed:
#
# - `fractiontopo` — fraction of subgrid sediment in each wet box,
# - `roughnesstopo` — proxy for topographic roughness.

# Loading the ETOPO data requires `Distances` and `NCDatasets` so the
# `AIBECSETOPOExt` extension activates.

using AIBECS, Plots
using Distances, NCDatasets
using JLD2 # required by `OCIM2.load`

# ## Load ETOPO and the target grid

grd, _ = OCIM2.load()
Z, lats, lons = ETOPO.load()                 # bedrock variant by default

# ## Fraction of sediment per box

f = fractiontopo(grd)
plothorizontalslice(f, grd, depth = 3000u"m", color = :viridis)

# Cells with `f > 0` at depth `z` correspond to grid boxes whose subgrid
# topography pokes above `z`.

# ## Roughness proxy

r = roughnesstopo(grd)
plothorizontalslice(r, grd, depth = 3000u"m", color = :inferno)

# Higher values flag boxes overlapping the seafloor with significant
# elevation variability — useful when parameterising boundary fluxes
# (e.g. sediment release) or interior mixing.
