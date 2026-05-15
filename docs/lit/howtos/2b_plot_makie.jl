#-------------------------------------------------------------------
# # [Plot basic things with Makie](@id plots-makie)
#-------------------------------------------------------------------

# This is the [Makie.jl](https://docs.makie.org/) companion to the
# [Plot basic things](@ref plots) guide. It shows the same examples,
# but built from the underlying data functions (`horizontalslice`,
# `meridionalslice`, `zonalmean`, `depthprofile` — re-exported from
# [OceanGrids](https://github.com/briochemc/OceanGrids.jl)) and
# rendered with Makie plotting calls instead of AIBECS' Plots.jl
# recipes.

# This guide is organized as follows
# - [Horizontal maps](@ref horizontal-plots-makie)
# - [Vertical slices](@ref vertical-plots-makie)
# - [Depth profiles](@ref profile-plots-makie)
# - [Makie extras](@ref makie-extras)

# We use the **CairoMakie** backend because the docs build runs
# headless in CI (no display server). For interactive work,
# `GLMakie` or `WGLMakie` use the same plotting API.

using AIBECS
using CairoMakie
using JLD2 # required by `OCIM2.load`
using Unitful
using Unitful: ustrip
CairoMakie.activate!(type = "png")

grd, _ = OCIM2.load()
dummy = cosd.(latvec(grd))

# Pre-strip the grid axis vectors once — Makie does not unwrap
# Unitful values automatically, so we hand it bare numbers and put
# the unit in the axis label.

lon = ustrip.(grd.lon)
lat = ustrip.(grd.lat)
depth = ustrip.(grd.depth)

#--------------------------------------------
# ## [Horizontal plots](@id horizontal-plots-makie)
#--------------------------------------------

# ### Horizontal slice

# To mirror `plothorizontalslice(dummy, grd, depth = 10)`, build the
# slice with `horizontalslice` and pass it to Makie's `heatmap`.
# Note the transpose: AIBECS' `horizontalslice` returns
# `(n_lat, n_lon)`, while Makie's `heatmap(x, y, z)` expects `z` to
# have shape `(length(x), length(y))`.

slice = AIBECS.horizontalslice(dummy, grd; depth = 10)
heatmap(lon, lat, slice';
    axis = (
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)",
        title = "dummy at 10 m",
    ),
)

# `horizontalslice` accepts the same `Unitful` depth values as the
# recipe-based version:

heatmap(lon, lat, AIBECS.horizontalslice(dummy, grd; depth = 10u"m")';
    axis = (xlabel = "Longitude (°E)", ylabel = "Latitude (°N)"))

# And different unit prefixes are converted under the hood by
# `horizontalslice` itself:

heatmap(lon, lat, AIBECS.horizontalslice(dummy, grd; depth = 3u"km")';
    axis = (xlabel = "Longitude (°E)", ylabel = "Latitude (°N)"))

# When the tracer carries units, we strip them at the call site and
# put the unit in the colorbar label. Makie's `rich` /
# `superscript` API gives typographic super/subscripts — nicer than
# relying on Unicode characters:

slice_units = AIBECS.horizontalslice(dummy * u"mol/m^3", grd; depth = 10u"m")
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Longitude (°E)", ylabel = "Latitude (°N)")
hm = heatmap!(ax, lon, lat, ustrip.(u"mol/m^3", slice_units)')
Colorbar(fig[1, 2], hm;
    label = rich("dummy (mol m", superscript("−3"), ")"))
fig

# The main trade-off vs. AIBECS' Plots.jl recipes is that Makie does
# not yet unwrap Unitful arrays automatically: we `ustrip` and label
# manually. That's a few extra characters per plot — but in exchange
# we get Makie's full layout system.

# Picking a colormap is a keyword on the plot call, just like with
# Plots.jl. The remaining axis/colorbar fine-tuning that
# `plot!(plt, xlabel = …, title = …, …)` does in the recipe example
# happens inline on the `Axis` and `Colorbar` constructors:

dummy .*= cosd.(lonvec(grd))
slice_balance = AIBECS.horizontalslice(dummy, grd; depth = 100)
fig = Figure()
ax = Axis(fig[1, 1];
    xlabel = "Lon",
    ylabel = "Lat",
    title = "The pacific as a whole",
)
hm = heatmap!(ax, lon, lat, slice_balance'; colormap = :balance)
Colorbar(fig[1, 2], hm; label = "dummy value")
fig

#----------------------------------------
# ## [Vertical plots](@id vertical-plots-makie)
#----------------------------------------

# ### Meridional slice

# `meridionalslice` returns `(n_depth, n_lat)`; we transpose for
# Makie and reverse the depth axis with `yreversed = true`:

dummy = cosd.(latvec(grd))
dummy .+= sqrt.(depthvec(grd)) / 30

mer = AIBECS.meridionalslice(dummy, grd; lon = 330)
heatmap(lat, depth, mer';
    axis = (
        xlabel = "Latitude (°N)",
        ylabel = "Depth (m)",
        yreversed = true,
        title = "Meridional slice at 330° E",
    ),
)

# ### Global zonal mean

# `zonalmean` already returns `(n_lat, n_depth)` (the AIBECS recipe
# transposes it before handing it to Plots), so no transpose is
# needed for Makie:

zmean = AIBECS.zonalmean(dummy, grd, 1)
heatmap(lat, depth, zmean;
    axis = (
        xlabel = "Latitude (°N)",
        ylabel = "Depth (m)",
        yreversed = true,
        title = "Global zonal mean",
    ),
)

# ### Basin-masked zonal mean

# As in the Plots.jl version, basin masks come from
# [OceanBasins.jl](https://github.com/briochemc/OceanBasins.jl):

using OceanBasins
OCEANS = oceanpolygons()
mPAC = ispacific(latvec(grd), lonvec(grd), OCEANS)
zmean_PAC = AIBECS.zonalmean(dummy, grd, mPAC)
heatmap(lat, depth, zmean_PAC;
    axis = (
        xlabel = "Latitude (°N)",
        ylabel = "Depth (m)",
        yreversed = true,
        title = "Pacific zonal mean",
    ),
)

# ### Another meridional slice

heatmap(lat, depth, AIBECS.meridionalslice(dummy, grd; lon = -30)';
    axis = (xlabel = "Latitude (°N)", ylabel = "Depth (m)",
            yreversed = true, title = "Meridional slice at -30° E"))

#----------------------------------------------------
# ## [Depth profiles](@id profile-plots-makie)
#----------------------------------------------------

# Depth profiles are 1-D — Makie's `lines` is the natural call. We
# put `depth` on the y-axis and flip it:

prof = AIBECS.depthprofile(dummy, grd; lonlat = (-30, 30))
lines(prof, depth;
    axis = (
        xlabel = "dummy",
        ylabel = "Depth (m)",
        yreversed = true,
        title = "Profile at (-30° E, 30° N)",
    ),
)

#----------------------------------------
# ## [Makie extras](@id makie-extras)
#----------------------------------------

# A few patterns Makie's layout system makes easy that are awkward
# with Plots.jl recipes.

# ### Two horizontal slices side by side with a shared colorbar

# Build a `Figure`, place two `Axis` objects in a 1×2 grid, plot
# into each, and share a single `Colorbar`:

dummy_h = cosd.(latvec(grd)) .* cosd.(lonvec(grd))
slice_surf = AIBECS.horizontalslice(dummy_h, grd; depth = 10u"m")
slice_deep = AIBECS.horizontalslice(dummy_h, grd; depth = 1000u"m")

cvals = filter(!isnan, vcat(vec(slice_surf), vec(slice_deep)))
vmin, vmax = extrema(cvals)

fig = Figure(size = (900, 350))
ax1 = Axis(fig[1, 1]; title = "10 m", xlabel = "Lon", ylabel = "Lat")
ax2 = Axis(fig[1, 2]; title = "1000 m", xlabel = "Lon")
hideydecorations!(ax2; grid = false)
hm1 = heatmap!(ax1, lon, lat, slice_surf'; colorrange = (vmin, vmax))
heatmap!(ax2, lon, lat, slice_deep'; colorrange = (vmin, vmax))
Colorbar(fig[1, 3], hm1; label = "dummy value")
fig

# ### Overlay a station marker on a horizontal slice

# Compositing with `scatter!` on the same `Axis` is a one-liner:

station_lon, station_lat = -30.0, 30.0
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Lon", ylabel = "Lat")
hm = heatmap!(ax, lon, lat, slice_surf')
scatter!(ax, [mod(station_lon, 360)], [station_lat];
    color = :red, markersize = 16, marker = :star5)
Colorbar(fig[1, 2], hm; label = "dummy value")
fig
