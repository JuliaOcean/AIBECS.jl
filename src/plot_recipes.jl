


"""
    HorizontalSlice(x, grd; depth)

Plots a horizontal slice of tracer `x` at depth `depth`.
"""
@userplot HorizontalSlice
@recipe function f(p::HorizontalSlice; depth=nothing)
    x, grd = p.args
    isnothing(depth) && error("you must specify the depth, e.g., `depth=100`")
    depth = convertdepth(depth)
    iz = finddepthindex(depth, grd)
    x3D = rearrange_into_3Darray(x, grd)
    @series begin
        seriestype := :contourf
        xguide --> "Longitude"
        yguide --> "Latitude"
        grd.lon, grd.lat, view(x3D, :, :, iz)
    end
end

#= TODO implement lon window for horizontal maps
lonshift(lon, lonstart) = mod(lon-lonstart, 360) + lonstart
lon0360(lon::Real) = mod(lon, 360)
lon0360(lon::Quantity) = mod(lon, 360u"°")
lon180W180E(lon::Real) = mod(lon + 180, 360) - 180
lon180W180E(lon::Quantity) = mod(lon + 180u"°", 360u"°") - 180u"°"
=#

"""
    finddepthindex(d, grd)

finds the depth index of `grd` closest to `depth`.
"""
finddepthindex(depth::Quantity, grd) = findmin(abs.(grd.depth .- depth))[2]
convertdepth(depth::Real) = depth * u"m"
convertdepth(depth::Quantity{U, Unitful.𝐋, V}) where {U,V} = depth
convertdepth(x) = error("Not a valid depth")

"""
    SurfaceMap(x, grd)

Plots a surface map of tracer `x`.
"""
@userplot SurfaceMap
@recipe function f(p::SurfaceMap)
    x, grd = p.args
    x3D = rearrange_into_3Darray(x, grd)
    lon, lat = grd.lon, grd.lat
    @series begin
        seriestype --> :contourf
        xguide --> "Longitude"
        yguide --> "Latitude"
        lon, lat, view(x3D, :, :, 1)
    end
end

"""
    VerticalIntegral(x, grd [; depthlim])

Plots the vertical integral of tracer `x`.
"""
@userplot VerticalIntegral
@recipe function f(p::VerticalIntegral; depthlim=extrema(p.args[2].depth))
    x, grd = p.args
    depthlim = convertdepth(depthlim)
    iz1, iz2 = finddepthindex(depthlim, grd)
    x3D = rearrange_into_3Darray(x, grd)[:,:,iz1:iz2]
    δz_3D = view(grd.δz_3D, :, :, iz1:iz2)
    x2D = nansum(x3D .* δz_3D, dims=3) ./ grd.wet3D[:,:,iz1]
    @series begin
        seriestype := :contourf
        xguide --> "Longitude"
        yguide --> "Latitude"
        grd.lon, grd.lat, view(x2D, :, :, 1)
    end
end
finddepthindex(depths::Tuple, grd) = ((finddepthindex(d, grd) for d in depths)...,)
convertdepth(depths::Tuple) = ((convertdepth(d) for d in depths)...,)

"""
    VerticalAverage(x, grd [; depthlim])

Plots the vertical average of tracer `x`.
"""
@userplot VerticalAverage
@recipe function f(p::VerticalAverage; depthlim=extrema(p.args[2].depth))
    x, grd = p.args
    depthlim = convertdepth(depthlim)
    iz1, iz2 = finddepthindex(depthlim, grd)
    x3D = rearrange_into_3Darray(x, grd)[:,:,iz1:iz2]
    v = vector_of_volumes(grd)
    v3D = rearrange_into_3Darray(v, grd)[:,:,iz1:iz2]
    x2D = nansum(x3D .* v3D, dims=3) ./ nansum(v3D, dims=3)
    @series begin
        seriestype := :contourf
        xguide --> "Longitude"
        yguide --> "Latitude"
        grd.lon, grd.lat, view(x2D, :, :, 1)
    end
end

# sun skipping nans
nansum(x; kwargs...) = sum(x .* (!isnan).(x); kwargs...)

"""
    ZonalSlice(x, grd; lon)

Plots a zonal slice of tracer `x` at longitude `lon`.
"""
@userplot ZonalSlice
@recipe function f(p::ZonalSlice; lon=nothing)
    isnothing(lon) && error("you must specify the longitude, e.g., `lon=100`")
    x, grd = p.args
    lon = convertlon(lon)
    ix = findlonindex(lon, grd)
    x3D = rearrange_into_3Darray(x, grd)
    @series begin
        seriestype := :contourf
        yflip := true
        ylims --> (0u"m", maximum(grd.depth))
        xguide --> "Latitude"
        yguide --> "Depth"
        grd.lat, grd.depth, permutedims(view(x3D,:,ix,:), [2,1])
    end
end
findlonindex(lon::Quantity, grd) = findmin(abs.(grd.lon .- lon))[2]
findlatindex(lat::Quantity, grd) = findmin(abs.(grd.lat .- lat))[2]
convertlon(x::Real) = mod(x * u"°", 360u"°")
convertlon(x::Quantity) = unit(x) == Unitful.° ? mod(x, 360u"°") : error("Not a valid lon")
convertlat(y::Real) = y * u"°"
convertlat(y::Quantity) = unit(y) == Unitful.° ? y : error("Not a valid lat")


"""
    MeridionalSlice(x, grd, lat)

Plots a Meridional slice of tracer `x` at longitude `lat`.
"""
@userplot MeridionalSlice
@recipe function f(p::MeridionalSlice; lat=nothing)
    isnothing(lat) && error("you must specify the latitude, e.g., `lat=-30`")
    x, grd = p.args
    lat = convertlat(lat)
    iy = findlatindex(lat, grd)
    x3D = rearrange_into_3Darray(x, grd)
    @series begin
        seriestype := :contourf
        yflip := true
        ylims --> (0u"m", maximum(grd.depth))
        xguide --> "Longitude"
        yguide --> "Depth"
        grd.lon, grd.depth, permutedims(view(x3D,iy,:,:), [2,1])
    end
end


"""
    ZonalAverage(x, grd; mask=1)

Plots a zonal average of tracer `x`.
"""
@userplot ZonalAverage
@recipe function f(p::ZonalAverage; mask=1)
    x, grd = p.args
    x3D = rearrange_into_3Darray(x .* mask, grd)
    v = vector_of_volumes(grd)
    v3D = rearrange_into_3Darray(v .* mask, grd)
    x2D = nansum(x3D .* v3D, dims=2) ./ nansum(v3D, dims=2)
    @series begin
        seriestype := :contourf
        yflip := true
        ylims --> (0u"m", maximum(grd.depth))
        xguide --> "Latitude"
        yguide --> "Depth"
        grd.lat, grd.depth, permutedims(view(x2D, :, 1, :), [2,1])
    end
end

"""
DepthProfile(x, grd; lonlat)

    Plots the profile of tracer `x` interpolated at `lonlat=(x,y)` coordinates.
"""
@userplot DepthProfile
@recipe function f(p::DepthProfile; lonlat=nothing)
    x, grd = p.args
    isnothing(lonlat) && error("you must specify the coordinates, e.g., `lonlat=(-30,30)`")
    lon, lat = lonlat
    lon, lat = convertlon(lon), convertlat(lat)
    x3D = rearrange_into_3Darray(x, grd)
    ndepth = length(grd.depth)
    knots = (grd.lat, grd.lon, 1:ndepth)
    itp = interpolate(knots, x3D, (Gridded(Linear()), Gridded(Linear()), NoInterp()))
    @series begin
        yflip --> true
        ylims --> (0u"m", maximum(grd.depth))
        yguide --> "Depth"
        itp(lat, lon, 1:ndepth), grd.depth
    end
end

# TODO fix this plot to smarlty reorder contour.
@userplot TransectContourfByDistance
@recipe function f(p::TransectContourfByDistance)
    x, grd, ct = p.args
    x3D = rearrange_into_3Darray(x, grd)
    depths = ustrip.(grd.depth)
    ndepths = length(depths)

    itp = interpolate(periodic_longitude(x3D), (BSpline(Linear()), BSpline(Linear()), NoInterp())) # interpolate linearly between the data points
    lats = range(ustrip(grd.lat[1]), ustrip(grd.lat[end]), length=length(grd.lat))
    lons = range(ustrip(grd.lon[1]), ustrip(grd.lon[1])+360, length=length(grd.lon)+1)
    stp = Interpolations.scale(itp, lats, lons, 1:ndepths) # re-scale to the actual domain
    etp = extrapolate(stp, (Line(), Periodic(), Line())) # periodic longitude

    n = length(ct)
    distances = cumsum(Distances.colwise(Haversine(6371.0),
                [view(ct.lon, :)'; view(ct.lat, :)'],
                [view(ct.lon, [1;1:n-1])'; view(ct.lat, [1;1:n-1])']))
    idx = unique(i -> distances[i], 1:length(ct)) # Remove stations at same (lat,lon)

    @series begin
        seriestype := :contourf
        yflip := true
        yticks --> Int.(round.(depths))
        ylims --> (0, maximum(depths))
        xguide --> "$(ct.name) distance (km)"
        yguide --> "Depth (m)"
        distances[idx], depths, [etp(lat, lon, i) for i in 1:ndepths, (lat, lon) in zip(ct.lat[idx], ct.lon[idx])]
    end
end

periodic_longitude(x3D::Array{T,3}) where T = view(x3D,:,[1:size(x3D,2); 1],:)
periodic_longitude(x2D::Array{T,2}) where T = view(x2D,:,[1:size(x3D,2); 1])


"""
    PlotCruiseTrack(ct; longitude_bounds)

Plots the cruise track `ct`.
"""
@userplot PlotCruiseTrack
@recipe function f(p::PlotCruiseTrack)
    ct = p.args[1]
    lons = [st.lon for st in ct.stations]
    lats = [st.lat for st in ct.stations]
    @series begin
        xguide --> "Longitude"
        yguide --> "Latitude"
        label --> ct.name
        markershape --> :hexagon
        linewidth --> 0
        linecolor --> :black
        markersize --> 3
        lons, lats
    end
end





"""
    MeridionalTransect(x, grd, ct)

Plots a Meridional transect of tracer `x` along cruise track `ct`.
"""
@userplot MeridionalTransect
@recipe function f(p::MeridionalTransect; ct=nothing)
    isnothing(ct) && error("you must include the cruise track with `ct=...`")
    x, grd = p.args
    x3D = rearrange_into_3Darray(x, grd)
    ndepth = length(grd.depth)

    itp = interpolate(periodic_longitude(x3D), (BSpline(Linear()), BSpline(Linear()), NoInterp())) # interpolate linearly between the data points

    lats = range(grd.lat[1], grd.lat[end], length=length(grd.lat))
    lons = range(grd.lon[1], grd.lon[1]+360u"°", length=length(grd.lon)+1)
    stp = Interpolations.scale(itp, lats, lons, 1:ndepth) # re-scale to the actual domain
    etp = extrapolate(stp, (Line(), Periodic(), Line())) # periodic longitude

    ctlats = [convertlat(st.lat) for st in ct.stations]
    ctlons = [convertlon(st.lon) for st in ct.stations]
    isort = sortperm(ctlats)
    idx = isort[unique(i -> ctlats[isort][i], 1:length(ct))] # Remove stations at same (lat,lon)
    @series begin
        seriestype := :contourf
        yflip := true
        ylims --> (0u"m", maximum(grd.depth))
        title --> "$(ct.name)"
        yguide --> "Depth"
        xguide --> "Latitude"
        ctlats[idx], grd.depth, [etp(lat, lon, i) for i in 1:ndepth, (lat, lon) in zip(ctlats[idx], ctlons[idx])]
    end
end


"""
    MeridionalScatterTransect(t)

Plots a scatter of the discrete obs of `t` in (lat,depth) space.
"""
@userplot MeridionalScatterTransect
@recipe function f(p::MeridionalScatterTransect)
    transect = p.args[1]
    depths = reduce(vcat, pro.depths for pro in transect.profiles)
    values = ustrip.(reduce(vcat, pro.values for pro in transect.profiles))
    lats = reduce(vcat, pro.station.lat * ones(length(pro)) for pro in transect.profiles)
    @series begin
        seriestype := :scatter
        yflip := true
        zcolor --> values
        markershape --> :circle
        xlim --> extrema(lats)
        clims --> (0, maximum(transect))
        label --> "$(transect.tracer) along $(transect.cruise)"
        yguide --> "Depth"
        xguide --> "Latitude"
        convertlat.(lats), convertdepth.(depths)
    end
end




"""
    ZonalTransect(x, grd, ct)

Plots a Zonal transect of tracer `x` along cruise track `ct`.
"""
@userplot ZonalTransect
@recipe function f(p::ZonalTransect; ct=nothing)
    isnothing(ct) && error("you must include the cruise track with `ct=...`")
    x, grd = p.args
    x3D = rearrange_into_3Darray(x, grd)
    ndepth = length(grd.depth)

    itp = interpolate(periodic_longitude(x3D), (BSpline(Linear()), BSpline(Linear()), NoInterp())) # interpolate linearly between the data points

    lats = range(grd.lat[1], grd.lat[end], length=length(grd.lat))
    lons = range(grd.lon[1], grd.lon[1]+360u"°", length=length(grd.lon)+1)
    stp = Interpolations.scale(itp, lats, lons, 1:ndepth) # re-scale to the actual domain
    etp = extrapolate(stp, (Line(), Periodic(), Line())) # periodic longitude

    ctlats = [convertlat(st.lat) for st in ct.stations]
    ctlons = [convertlon(st.lon) for st in ct.stations]
    isort = sortperm(ctlons)
    idx = isort[unique(i -> ctlons[isort][i], 1:length(ct))] # Remove stations at same (lat,lon)
    @series begin
        seriestype := :contourf
        yflip := true
        ylims --> (0u"m", maximum(grd.depth))
        title --> "$(ct.name)"
        yguide --> "Depth"
        xguide --> "Longitude"
        ctlons[idx], grd.depth, [etp(lat, lon, i) for i in 1:ndepth, (lat, lon) in zip(ctlats[idx], ctlons[idx])]
    end
end


"""
    ZonalScatterTransect(t)

Plots a scatter of the discrete obs of `t` in (lat,depth) space.
"""
@userplot ZonalScatterTransect
@recipe function f(p::ZonalScatterTransect)
    transect = p.args[1]
    depths = reduce(vcat, pro.depths for pro in transect.profiles)
    values = ustrip.(reduce(vcat, pro.values for pro in transect.profiles))
    lons = reduce(vcat, pro.station.lon * ones(length(pro)) for pro in transect.profiles)
    @series begin
        seriestype := :scatter
        zcolor := values
        yflip := true
        markershape --> :circle
        xlim --> extrema(lons)
        clims --> (0, maximum(transect))
        label --> "$(transect.tracer) along $(transect.cruise)"
        yguide --> "Depth"
        xguide --> "Longitude"
        colorbar_title --> string(unit(transect))
        convertlon.(lons), convertdepth.(depths)
    end
end





"""
    RatioAtStation(x, y, grd, station, depthlims=(0,Inf))

Plots a meridional transect of tracer `x` along cruise track `ct`.

The keyword argument `zlims=(ztop, zbottom)` can be provided
if you only want to only plot for depths `z ∈ (ztop, zbottom)`.
(`z` is positive downwards in this case)
"""
@userplot RatioAtStation
@recipe function f(p::RatioAtStation; depthlims=(0,Inf))
    depthlims = convertdepth.(depthlims)
    x, y, grd, st = p.args
    x3D = rearrange_into_3Darray(x, grd)
    y3D = rearrange_into_3Darray(y, grd)
    ndepth = length(grd.depth)
    knots = (grd.lat, grd.lon, 1:ndepth)
    itpx = interpolate(knots, x3D, (Gridded(Linear()), Gridded(Linear()), NoInterp()))
    itpy = interpolate(knots, y3D, (Gridded(Linear()), Gridded(Linear()), NoInterp()))
    iz = findall(depthlims[1] .≤ grd.depth .≤ depthlims[2])
    ibot = findfirst(grd.depth .> depthlims[2])
    x1D = itpx(convertlat(st.lat), mod(convertlon(st.lon), 360u"°"), iz)
    y1D = itpy(convertlat(st.lat), mod(convertlon(st.lon), 360u"°"), iz)
    @series begin
        seriestype := :scatter
        label --> st.name
        markershape --> :circle
        seriescolor --> :deep
        marker_z --> -ustrip(grd.depth) # TODO fix in UnitfulRecipes
        colorbar_title --> "Depth (m)" # TODO fix in UniftulRecipes
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




#=




=#