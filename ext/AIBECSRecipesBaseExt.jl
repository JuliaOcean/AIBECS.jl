module AIBECSRecipesBaseExt

using AIBECS
using RecipesBase
using Interpolations
using OceanGrids
using Unitful
using Unitful: °, ustrip

import AIBECS: ratioatstation, ratioatstation!, plotstencil, plotstencil!

const lonticks = (-180:30:360, ["180°", "", "120°W", "", "60°W", "", "0°", "", "60°E", "", "120°E", "", "180°", "", "120°W", "", "60°W", "", "0°"])
const latticks = (-90:30:90, ["SP", "60°S", "30°S", "EQ", "30°N", "60°N", "NP"])
function prettyticks(axis, defaultticks)
    idx = findall(ustrip(minimum(axis)) .≤ defaultticks[1] .≤ ustrip(maximum(axis)))
    return (defaultticks[1][idx], defaultticks[2][idx])
end

#============================
horizontal maps / (x,y) plots
============================#
@userplot HorizontalPlane
@recipe function f(p::HorizontalPlane)
    x, y, z = p.args
    @series begin
        seriestype --> :heatmap
        size --> (800, 400)
        aspect_ratio --> :equal
        framestyle --> :box
        xguide --> ""
        yguide --> ""
        xticks --> prettyticks(x, lonticks)
        yticks --> prettyticks(y, latticks)
        x, y, z
    end
end
AIBECS.plothorizontalslice(x, grd; depth = nothing, kwargs...) = horizontalplane(grd.lon, grd.lat, AIBECS.horizontalslice(x, grd; depth = depth); kwargs...)
AIBECS.plothorizontalslice!(x, grd; depth = nothing, kwargs...) = horizontalplane!(grd.lon, grd.lat, AIBECS.horizontalslice(x, grd; depth = depth); kwargs...)
AIBECS.plothorizontalslice!(plt, x, grd; depth = nothing, kwargs...) = horizontalplane!(plt, grd.lon, grd.lat, AIBECS.horizontalslice(x, grd; depth = depth); kwargs...)
AIBECS.surfacemap(args...; kwargs...) = AIBECS.plothorizontalslice(args...; kwargs..., depth = 0)
AIBECS.surfacemap!(args...; kwargs...) = AIBECS.plothorizontalslice!(args...; kwargs..., depth = 0)
AIBECS.plot∫dz(x, grd; mask = 1, kwargs...) = horizontalplane(grd.lon, grd.lat, AIBECS.∫dz(x, grd, mask); kwargs...)
AIBECS.plot∫dz!(plt, x, grd; mask = 1, kwargs...) = horizontalplane!(plt, grd.lon, grd.lat, AIBECS.∫dz(x, grd, mask); kwargs...)
AIBECS.plotverticalmean(x, grd; mask = 1, kwargs...) = horizontalplane(grd.lon, grd.lat, AIBECS.verticalmean(x, grd, mask); kwargs...)

@userplot MiniMap
@recipe function f(p::MiniMap; central_longitude = 200°)
    grd, = p.args
    wlon = central_longitude - 180°
    xgrd = @. mod(grd.lon - wlon, 360°) + wlon
    ix = sortperm(xgrd)
    @series begin
        seriestype --> :heatmap
        seriescolor --> :GnBu_3
        rev --> true
        xguide --> ""
        yguide --> ""
        framestyle  --> false
        showaxis --> false
        legend --> :none
        colorbar --> false
        grid --> false
        xgrd[ix], grd.lat, view(AIBECS.iswet(grd), :, ix, 1)
    end
end

#================================
Vertical–meridional / (y,z) plots
================================#
@userplot MeridionalPlane
@recipe function f(p::MeridionalPlane)
    y, z, v = p.args
    @series begin
        seriestype --> :heatmap
        framestyle --> :box
        xguide --> ""
        yguide --> "Depth"
        xticks --> prettyticks(y, latticks)
        yflip --> true
        y, z, v
    end
end
AIBECS.plotmeridionalslice(x, grd; lon = nothing, kwargs...) = meridionalplane(grd.lat, grd.depth, AIBECS.meridionalslice(x, grd; lon = lon); kwargs...)
AIBECS.plotmeridionalslice!(x, grd; lon = nothing, kwargs...) = meridionalplane!(grd.lat, grd.depth, AIBECS.meridionalslice(x, grd; lon = lon); kwargs...)
AIBECS.plotmeridionalslice!(plt, x, grd; lon = nothing, kwargs...) = meridionalplane!(plt, grd.lat, grd.depth, AIBECS.meridionalslice(x, grd; lon = lon); kwargs...)
AIBECS.plotzonalmean(x, grd; mask = 1, kwargs...) = meridionalplane(grd.lat, grd.depth, AIBECS.zonalmean(x, grd, mask)'; kwargs...)
AIBECS.plotzonalmean!(x, grd; mask = 1, kwargs...) = meridionalplane!(grd.lat, grd.depth, AIBECS.zonalmean(x, grd, mask)'; kwargs...)
AIBECS.plotzonalmean!(plt, x, grd; mask = 1, kwargs...) = meridionalplane!(plt, grd.lat, grd.depth, AIBECS.zonalmean(x, grd, mask)'; kwargs...)
AIBECS.plot∫dx(x, grd; mask = 1, kwargs...) = meridionalplane(grd.lat, grd.depth, AIBECS.∫dx(x, grd, mask)'; kwargs...)

#===========================
Vertical–zonal / (x,z) plots
===========================#
@userplot ZonalPlane
@recipe function f(p::ZonalPlane)
    x, z, v = p.args
    @series begin
        seriestype --> :heatmap
        xguide --> "Longitude"
        yguide --> "Depth"
        yflip --> true
        x, z, v
    end
end
AIBECS.plotzonalslice(x, grd; lat = nothing, kwargs...) = zonalplane(grd.lon, grd.depth, AIBECS.zonalslice(x, grd; lat = lat); kwargs...)
AIBECS.plotmeridionalmean(x, grd; mask = 1, kwargs...) = zonalplane(grd.lon, grd.depth, AIBECS.meridionalmean(x, grd, mask)'; kwargs...)
AIBECS.plot∫dy(x, grd; mask = 1, kwargs...) = zonalplane(grd.lon, grd.depth, AIBECS.∫dy(x, grd, mask)'; kwargs...)

#=================
Vertical / z plots
=================#
@userplot ZPlot
@recipe function f(p::ZPlot)
    x, z = p.args
    @series begin
        yguide --> "Depth"
        yflip --> true
        xguide_position --> :top
        xmirror --> :true
        x, z
    end
end
AIBECS.plot∫dxdy(x, grd; mask = 1, kwargs...) = zplot(AIBECS.∫dxdy(x, grd, mask), grd.depth; kwargs...)
AIBECS.plot∫dxdy!(x, grd; mask = 1, kwargs...) = zplot!(AIBECS.∫dxdy(x, grd, mask), grd.depth; kwargs...)
AIBECS.plot∫dxdy!(plt, x, grd; mask = 1, kwargs...) = zplot!(plt, AIBECS.∫dxdy(x, grd, mask), grd.depth; kwargs...)
AIBECS.plothorizontalmean(x, grd; mask = 1, kwargs...) = zplot(AIBECS.horizontalmean(x, grd, mask), grd.depth; kwargs...)
AIBECS.plothorizontalmean!(x, grd; mask = 1, kwargs...) = zplot!(AIBECS.horizontalmean(x, grd, mask), grd.depth; kwargs...)
AIBECS.plothorizontalmean!(plt, x, grd; mask = 1, kwargs...) = zplot!(plt, AIBECS.horizontalmean(x, grd, mask), grd.depth; kwargs...)
AIBECS.plotdepthprofile(x, grd; lonlat = nothing, kwargs...) = zplot(AIBECS.depthprofile(x, grd; lonlat = lonlat), grd.depth; kwargs...)
AIBECS.plotdepthprofile!(plt, x, grd; lonlat = nothing, kwargs...) = zplot!(plt, AIBECS.depthprofile(x, grd; lonlat = lonlat), grd.depth; kwargs...)

#======================================================
Vertical / (dist,z) plots (for transect/section slices)
======================================================#
@userplot VerticalPlane
@recipe function f(p::VerticalPlane)
    dist, z, v = p.args
    @series begin
        seriestype --> :heatmap
        xguide --> "Distance"
        yguide --> "Depth"
        yflip --> true
        grid --> false
        dist, z, v
    end
end
AIBECS.plottransect(x, grd; ct = nothing, start = :south, show_stations = false, insetmap = false, kwargs...) = verticalplane(AIBECS.verticalsection2(x, grd; ct = ct)...; kwargs...)

@userplot RatioAtStation
@recipe function f(p::RatioAtStation; depthlims = (0, Inf))
    depthlims = OceanGrids.convertdepth.(depthlims)
    x, y, grd, st = p.args
    x3D = AIBECS.rearrange_into_3Darray(x, grd)
    y3D = AIBECS.rearrange_into_3Darray(y, grd)
    ndepth = length(grd.depth)
    knots = (grd.lat, grd.lon, 1:ndepth)
    itpx = interpolate(knots, x3D, (Gridded(Linear()), Gridded(Linear()), NoInterp()))
    itpy = interpolate(knots, y3D, (Gridded(Linear()), Gridded(Linear()), NoInterp()))
    iz = findall(depthlims[1] .≤ grd.depth .≤ depthlims[2])
    ibot = findfirst(grd.depth .> depthlims[2])
    x1D = itpx(AIBECS.convertlat(st.lat), mod(AIBECS.convertlon(st.lon), 360°), iz)
    y1D = itpy(AIBECS.convertlat(st.lat), mod(AIBECS.convertlon(st.lon), 360°), iz)
    @series begin
        seriestype := :scatter
        label --> st.name
        markershape --> :circle
        seriescolor --> :deep
        marker_z --> ustrip.(-grd.depth)
        colorbar_title --> "Depth (m)"
        markersize --> 4
        legend --> :topleft
        x1D, y1D
    end
end

#======================================================
Diagnostic recipes
======================================================#
@userplot PlotStencil
@recipe function f(p::PlotStencil)
    st, = p.args
    @series begin
        linecolor --> :black
        label --> ""
        convert.(Tuple, vcat([[0 * x, x, 0 * x] for x in st]...))
    end
    @series begin
        marker_z --> [x.I[3] for x in st]
        label --> ""
        xguide --> "lat index"
        yguide --> "lon index"
        zguide --> "depth index"
        seriescolor --> :darkrainbow
        linewidth --> 0
        markershape --> :circle
        convert.(Tuple, st)
    end
end

end # module
