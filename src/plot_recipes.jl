

#============================
horizontal maps / (x,y) plots
============================#
@userplot HorizontalPlane
@recipe function f(p::HorizontalPlane)
    x, y, z = p.args
    @series begin
        seriestype --> :heatmap
        xguide --> "Longitude"
        yguide --> "Latitude"
        x, y, z
    end
end
"""
    plothorizontalslice(x, grd; depth)

Plots a heatmap of the horizontal slice of tracer `x` at depth `depth`.
"""
plothorizontalslice(x, grd; depth=nothing, kwargs...) = horizontalplane(grd.lon, grd.lat, horizontalslice(x, grd; depth=depth); kwargs...)
plothorizontalslice!(x, grd; depth=nothing, kwargs...) = horizontalplane!(grd.lon, grd.lat, horizontalslice(x, grd; depth=depth); kwargs...)
plothorizontalslice!(plt, x, grd; depth=nothing, kwargs...) = horizontalplane!(plt, grd.lon, grd.lat, horizontalslice(x, grd; depth=depth); kwargs...)
"""
    surfacemap(x, grd)

Plots a surface heatmap of tracer `x`.
"""
surfacemap(args...; kwargs...) = plothorizontalslice(args...; kwargs..., depth=0)
surfacemap!(args...; kwargs...) = plothorizontalslice!(args...; kwargs..., depth=0)
"""
    plotverticalintegral(x, grd, mask=1)

Plots the vertical integral of tracer `x`.
"""
plot∫dz(x, grd; mask=1, kwargs...) = horizontalplane(grd.lon, grd.lat, ∫dz(x, grd, mask); kwargs...)
plot∫dz!(plt, x, grd; mask=1, kwargs...) = horizontalplane!(plt, grd.lon, grd.lat, ∫dz(x, grd, mask); kwargs...)
plotverticalintegral = plot∫dz
plotverticalintegral! = plot∫dz!
"""
    plotverticalaverage(x, grd, mask=1)

Plots the vertical average of tracer `x`.
"""
plotverticalmean(x, grd; mask=1, kwargs...) = horizontalplane(grd.lon, grd.lat, verticalmean(x, grd, mask); kwargs...)
plotverticalaverage = plotverticalmean
export plotverticalaverage, plotverticalmean, plotverticalintegral, plotverticalintegral!, plot∫dz, plot∫dz!, surfacemap, surfacemap!, plothorizontalslice

"""
    minimap(grd; central_longitude=200°)

Plots a surface map of grd with no ticks, labels, or colorbar.

The goal of this function is to provide a way to plot both a minimap and a
cruise track when plotting a transect:
```
plottransect(dummy, grd; ct=sort(ct))
plot!(inset_subplots=bbox(0.73,0.73,0.15,0.15), subplot=2)
minimap!(grd; subplot=2)
cruisetrackplots!(sort(ct), subplot=2)
```
"""
@userplot MiniMap
@recipe function f(p::MiniMap; central_longitude=200°)
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
        xgrd[ix], grd.lat, view(iswet(grd),:,ix,1)
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
        xguide --> "Latitude"
        yguide --> "Depth"
        yflip --> true
        y, z, v
    end
end
"""
    plotmeridionalslice(x, grd; lon)

Plots a meridional slice of tracer `x` at longitude `lon`.
"""
plotmeridionalslice(x, grd; lon=nothing, kwargs...) = meridionalplane(grd.lat, grd.depth, meridionalslice(x, grd; lon=lon); kwargs...)
"""
    plotzonalmean(x, grd; mask=1)

Plots a zonal average of tracer `x`.
"""
plotzonalmean(x, grd; mask=1, kwargs...) = meridionalplane(grd.lat, grd.depth, zonalmean(x, grd, mask)'; kwargs...)
plotzonalaverage = plotzonalmean
"""
    plotzonalintegral(x, grd; mask=1)

Plots a zonal integral of tracer `x`.
"""
plot∫dx(x, grd; mask=1, kwargs...) = meridionalplane(grd.lat, grd.depth, ∫dx(x, grd, mask)'; kwargs...)
plotzonalintegral = plot∫dx
export plotzonalintegral, plot∫dx, plotzonalmean, plotzonalaverage, plotmeridionalslice



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
"""
    plotzonalslice(x, grd; lat)

Plots a zonal slice of tracer `x` at latitude `lat`.
"""
plotzonalslice(x, grd; lat=nothing, kwargs...) = zonalplane(grd.lon, grd.depth, zonalslice(x, grd; lat=lat); kwargs...)
"""
    plotmeridionalmean(x, grd; mask=1)

Plots a meridional average of tracer `x`.
"""
plotmeridionalmean(x, grd; mask=1, kwargs...) = zonalplane(grd.lon, grd.depth, meridionalmean(x, grd, mask)'; kwargs...)
plotmeridionalaverage = plotmeridionalmean
"""
    plotmeridionalintegral(x, grd; mask=1)

Plots a meridional integral of tracer `x`.
"""
plot∫dy(x, grd; mask=1, kwargs...) = zonalplane(grd.lon, grd.depth, ∫dy(x, grd, mask)'; kwargs...)
plotmeridionalintegral = plot∫dy
export plotmeridionalintegral, plot∫dy, plotmeridionalmean, plotmeridionalaverage, plotzonalslice

#=================
Vertical / z plots
=================#
@userplot ZPlot
@recipe function f(p::ZPlot)
    x, z = p.args
    @series begin
        yguide --> "Depth"
        yflip --> true
        x, z
    end
end

"""
    horizontalintegral(x, grd; mask=1)

Plots a horizontal integral of tracer `x`.
"""
plot∫dxdy(x, grd; mask=1, kwargs...) = zplot(∫dxdy(x, grd, mask), grd.depth; kwargs...)
plot∫dxdy!(x, grd; mask=1, kwargs...) = zplot!(∫dxdy(x, grd, mask), grd.depth; kwargs...)
plot∫dxdy!(plt, x, grd; mask=1, kwargs...) = zplot!(plt, ∫dxdy(x, grd, mask), grd.depth; kwargs...)
plothorizontalintegral = plot∫dxdy
plothorizontalintegral! = plot∫dxdy!
"""
    plothorizontalaverage(x, grd; mask=1)

Plots a horizontal average of tracer `x`.
"""
plothorizontalmean(x, grd; mask=1, kwargs...) = zplot(horizontalmean(x, grd, mask), grd.depth; kwargs...)
plothorizontalmean!(plt, x, grd; mask=1, kwargs...) = zplot!(plt, horizontalmean(x, grd, mask), grd.depth; kwargs...)
plothorizontalaverage = plothorizontalmean
plothorizontalaverage! = plothorizontalmean!
"""
    plotdepthprofile(x, grd; lonlat)

Plots the profile of tracer `x` interpolated at `lonlat=(x,y)` coordinates.
"""
plotdepthprofile(x, grd; lonlat=nothing, kwargs...) = zplot(depthprofile(x, grd; lonlat=lonlat), grd.depth; kwargs...)
export plot∫dxdy, plot∫dxdy!, plothorizontalintegral, plothorizontalintegral!, plothorizontalmean, plothorizontalmean!, plothorizontalaverage, plothorizontalaverage!, plotdepthprofile

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
"""
    plottransect(x, grd; ct=ct)

Plots the transect of tracer `x` along transect `ct`.
"""
plottransect(x, grd; ct=nothing, start=:south, show_stations=false, insetmap=false, kwargs...) = verticalplane(verticalsection(x, grd; ct=ct)...; kwargs...)
export plottransect




"""
    RatioAtStation(x, y, grd, station, depthlims=(0,Inf))

Plots a meridional transect of tracer `x` along cruise track `ct`.

The keyword argument `zlims=(ztop, zbottom)` can be provided
if you only want to only plot for depths `z ∈ (ztop, zbottom)`.
(`z` is positive downwards in this case)
"""
@userplot RatioAtStation
@recipe function f(p::RatioAtStation; depthlims=(0,Inf))
    depthlims = OceanGrids.convertdepth.(depthlims)
    x, y, grd, st = p.args
    x3D = rearrange_into_3Darray(x, grd)
    y3D = rearrange_into_3Darray(y, grd)
    ndepth = length(grd.depth)
    knots = (grd.lat, grd.lon, 1:ndepth)
    itpx = interpolate(knots, x3D, (Gridded(Linear()), Gridded(Linear()), NoInterp()))
    itpy = interpolate(knots, y3D, (Gridded(Linear()), Gridded(Linear()), NoInterp()))
    iz = findall(depthlims[1] .≤ grd.depth .≤ depthlims[2])
    ibot = findfirst(grd.depth .> depthlims[2])
    x1D = itpx(convertlat(st.lat), mod(convertlon(st.lon), 360°), iz)
    y1D = itpy(convertlat(st.lat), mod(convertlon(st.lon), 360°), iz)
    @series begin
        seriestype := :scatter
        label --> st.name
        markershape --> :circle
        seriescolor --> :deep
        marker_z --> ustrip.(-grd.depth)
        colorbar_title --> "Depth (m)" # TODO fix in UnitfulRecipes
        markersize --> 4
        legend --> :topleft
        x1D, y1D
    end
end



"""
    PlotParameter(p::AbstractParameters, s)

Plots the PDF of parameter `p` with symbol `s`
"""
@userplot PlotParameter
@recipe function f(plt::PlotParameter)
    d, v, initv, s, u = extract_dvisu(plt.args...)
    xs = default_range(d)
    ys = [pdf(d, x) for x in xs] # the PDF
    xu = xs * upreferred(u) .|> u .|> ustrip
    vu = v * upreferred(u) |> u |> ustrip
    @series begin
        label --> "prior"
        xu, ys
    end
    @series begin
        label --> "initial value"
        initv * [1,1], [0, pdf(d, v)]
    end
    @series begin
        xguide --> "$s ($(string(u)))"
        yguide --> "PDF"
        label --> "value"
        markershape --> :o
        [vu], [pdf(d, v)]
    end
end

"""
    PlotParameters(p::AbstractParameters)

Plots the PDF of all the flattenable parameters in `p`.
"""
@userplot PlotParameters
@recipe function f(plt::PlotParameters)
    p, = plt.args
    layout --> (length(p),1)
    for (i,s) in enumerate(flattenable_symbols(p))
        d, v, initvu, s, u = extract_dvisu(p,s)
        xs = default_range(d)
        ys = [pdf(d, x) for x in xs] # the PDF
        xu = xs * upreferred(u) .|> u .|> ustrip
        vu = v * upreferred(u) |> u |> ustrip
        initv = ustrip(upreferred(initvu * units(p,s)))
        @series begin
            label --> "prior"
            subplot := i
            xu, ys
        end
        @series begin
            label --> "initial value"
            subplot := i
            initvu * [1,1], [0, pdf(d, initv)]
        end
        @series begin
            xguide --> "$s ($(string(u)))"
            yguide --> "PDF"
            label --> "value"
            markershape --> :o
            subplot := i
            [vu], [pdf(d, v)]
        end
    end
end

extract_dvisu(d, v, initv, s, u) = d, v, initv, s, u
extract_dvisu(p::AbstractParameters, s) = prior(p,s), UnPack.unpack(p, Val(s)), initial_value(p,s), s, units(p,s)


function default_range(d::Distribution, n = 4)
    μ, σ = mean(d), std(d)
    xmin = max(μ - n*σ, minimum(support(d)))
    xmax = min(μ + n*σ, maximum(support(d)))
    range(xmin, xmax, length=1001)
end



#======================================================
Diagnostic recipes
======================================================#
"""
    PlotStencil(st)

Plots the stencil `st`.
"""
@userplot PlotStencil
@recipe function f(p::PlotStencil)
    st, = p.args
    @series begin
        linecolor --> :black
        label --> ""
        convert.(Tuple, vcat([[0*x, x, 0*x] for x in st]...))
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

