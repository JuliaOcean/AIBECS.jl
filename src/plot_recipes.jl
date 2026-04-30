# Plotting recipes are provided by the AIBECSRecipesBaseExt and
# AIBECSRecipesParametersExt extensions. The functions below are stubs;
# methods are added when the user loads the trigger packages
# (e.g. `using Plots` for the spatial recipes, plus `using Distributions`
# for parameter PDF recipes).

"""
    plothorizontalslice(x, grd; depth, kwargs...)

Plot a horizontal (lon-lat) slice of tracer `x` at the given `depth`.
Requires `using Plots`.
"""
function plothorizontalslice end

"""
    plothorizontalslice!(plt, x, grd; depth, kwargs...)

In-place version of [`plothorizontalslice`](@ref).
Requires `using Plots`.
"""
function plothorizontalslice! end

"""
    surfacemap(x, grd; kwargs...)

Plot a horizontal map of tracer `x` at the surface. Requires `using Plots`.
"""
function surfacemap end

"""
    surfacemap!(plt, x, grd; kwargs...)

In-place version of [`surfacemap`](@ref). Requires `using Plots`.
"""
function surfacemap! end

"""
    plot∫dz(x, grd; kwargs...)

Plot the depth-integral of tracer `x` as a horizontal map.
Aliased as `plotverticalintegral`. Requires `using Plots`.
"""
function plot∫dz end

"""
    plot∫dz!(plt, x, grd; kwargs...)

In-place version of [`plot∫dz`](@ref). Requires `using Plots`.
"""
function plot∫dz! end
const plotverticalintegral = plot∫dz
const plotverticalintegral! = plot∫dz!

"""
    plotverticalmean(x, grd; kwargs...)

Plot the depth-weighted mean of tracer `x` as a horizontal map.
Aliased as `plotverticalaverage`. Requires `using Plots`.
"""
function plotverticalmean end
const plotverticalaverage = plotverticalmean

"""
    plotmeridionalslice(x, grd; lon, kwargs...)

Plot a meridional (lat-depth) slice of tracer `x` at longitude `lon`.
Requires `using Plots`.
"""
function plotmeridionalslice end

"""
    plotmeridionalslice!(plt, x, grd; lon, kwargs...)

In-place version of [`plotmeridionalslice`](@ref). Requires `using Plots`.
"""
function plotmeridionalslice! end

"""
    plotzonalmean(x, grd; kwargs...)

Plot the zonally averaged tracer `x` as a lat-depth map.
Aliased as `plotzonalaverage`. Requires `using Plots`.
"""
function plotzonalmean end

"""
    plotzonalmean!(plt, x, grd; kwargs...)

In-place version of [`plotzonalmean`](@ref). Requires `using Plots`.
"""
function plotzonalmean! end
const plotzonalaverage = plotzonalmean
const plotzonalaverage! = plotzonalmean!

"""
    plot∫dx(x, grd; kwargs...)

Plot the zonal integral of tracer `x` as a lat-depth map.
Aliased as `plotzonalintegral`. Requires `using Plots`.
"""
function plot∫dx end
const plotzonalintegral = plot∫dx

"""
    plotzonalslice(x, grd; lat, kwargs...)

Plot a zonal (lon-depth) slice of tracer `x` at latitude `lat`.
Requires `using Plots`.
"""
function plotzonalslice end

"""
    plotmeridionalmean(x, grd; kwargs...)

Plot the meridionally averaged tracer `x` as a lon-depth map.
Aliased as `plotmeridionalaverage`. Requires `using Plots`.
"""
function plotmeridionalmean end
const plotmeridionalaverage = plotmeridionalmean

"""
    plot∫dy(x, grd; kwargs...)

Plot the meridional integral of tracer `x` as a lon-depth map.
Aliased as `plotmeridionalintegral`. Requires `using Plots`.
"""
function plot∫dy end
const plotmeridionalintegral = plot∫dy

"""
    plot∫dxdy(x, grd; kwargs...)

Plot the horizontally integrated tracer `x` as a depth profile.
Aliased as `plothorizontalintegral`. Requires `using Plots`.
"""
function plot∫dxdy end

"""
    plot∫dxdy!(plt, x, grd; kwargs...)

In-place version of [`plot∫dxdy`](@ref). Requires `using Plots`.
"""
function plot∫dxdy! end
const plothorizontalintegral = plot∫dxdy
const plothorizontalintegral! = plot∫dxdy!

"""
    plothorizontalmean(x, grd; kwargs...)

Plot the area-weighted horizontal mean of tracer `x` as a depth profile.
Aliased as `plothorizontalaverage`. Requires `using Plots`.
"""
function plothorizontalmean end

"""
    plothorizontalmean!(plt, x, grd; kwargs...)

In-place version of [`plothorizontalmean`](@ref). Requires `using Plots`.
"""
function plothorizontalmean! end
const plothorizontalaverage = plothorizontalmean
const plothorizontalaverage! = plothorizontalmean!

"""
    plotdepthprofile(x, grd; lon, lat, kwargs...)

Plot the depth profile of tracer `x` at the given `lon` and `lat`.
Requires `using Plots`.
"""
function plotdepthprofile end

"""
    plotdepthprofile!(plt, x, grd; lon, lat, kwargs...)

In-place version of [`plotdepthprofile`](@ref). Requires `using Plots`.
"""
function plotdepthprofile! end

"""
    plottransect(x, grd, ct; kwargs...)

Plot tracer `x` along the geographic transect `ct`. Requires `using Plots`.
"""
function plottransect end

"""
    ratioatstation(x, y, grd, station; kwargs...)

Scatter plot of `x` against `y` for the column at `station`,
useful for tracer-vs-tracer property plots. Requires `using Plots`.
"""
function ratioatstation end

"""
    ratioatstation!(plt, x, y, grd, station; kwargs...)

In-place version of [`ratioatstation`](@ref). Requires `using Plots`.
"""
function ratioatstation! end

"""
    plotstencil(T, grd; kwargs...)

Visualize the [`stencil`](@ref) of transport operator `T` on grid `grd`.
Requires `using Plots`.
"""
function plotstencil end

"""
    plotstencil!(plt, T, grd; kwargs...)

In-place version of [`plotstencil`](@ref). Requires `using Plots`.
"""
function plotstencil! end

export plothorizontalslice, plothorizontalslice!, surfacemap, surfacemap!,
    plot∫dz, plot∫dz!, plotverticalintegral, plotverticalintegral!,
    plotverticalmean, plotverticalaverage,
    plotmeridionalslice, plotmeridionalslice!,
    plotzonalmean, plotzonalmean!, plotzonalaverage, plotzonalaverage!,
    plot∫dx, plotzonalintegral,
    plotzonalslice, plotmeridionalmean, plotmeridionalaverage,
    plot∫dy, plotmeridionalintegral,
    plot∫dxdy, plot∫dxdy!, plothorizontalintegral, plothorizontalintegral!,
    plothorizontalmean, plothorizontalmean!, plothorizontalaverage, plothorizontalaverage!,
    plotdepthprofile, plotdepthprofile!,
    plottransect,
    ratioatstation, ratioatstation!,
    plotstencil, plotstencil!
