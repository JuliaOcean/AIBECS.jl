"""
    HorizontalSlice(x, grd, depth)

Plots a horizontal slice of tracer `x` at depth `depth`.
"""
@userplot HorizontalSlice
@recipe function f(p::HorizontalSlice)
    x, grd, depth = p.args
    xunit = string(unit(x[1]))
    x3D = rearrange_into_3Darray(ustrip.(x), grd)
    lon, lat = grd.lon .|> ustrip, grd.lat .|> ustrip
    iz = findfirst(ustrip.(grd.depth) .≥ ustrip(upreferred(depth)))
    isnothing(iz) && (iz = length(grd.depth))
    @series begin
        seriestype := :contourf
        xlabel := "Longitude"
        ylabel := "Latitude"
        colorbar_title := xunit
        lon, lat, view(x3D, :, :, iz)
    end
end

"""
    SurfaceMap(x, grd)

Plots a surface map of tracer `x`.
"""
@userplot SurfaceMap
@recipe function f(p::SurfaceMap)
    x, grd = p.args
    xunit = string(unit(x[1]))
    x3D = rearrange_into_3Darray(ustrip.(x), grd)
    lon, lat = grd.lon .|> ustrip, grd.lat .|> ustrip
    @series begin
        seriestype := :contourf
        xlabel := "Longitude"
        ylabel := "Latitude"
        colorbar_title := xunit
        lon, lat, view(x3D, :, :, 1)
    end
end

"""
    VerticalIntegral(x, grd)

Plots the vertical integral of tracer `x`.
"""
@userplot VerticalIntegral
@recipe function f(p::VerticalIntegral)
    x, grd = p.args
    intunit = string(unit(x[1]) * u"m")
    x3D = rearrange_into_3Darray(ustrip.(x), grd)
    lon, lat = grd.lon .|> ustrip, grd.lat .|> ustrip
    δz_3D = ustrip.(grd.δz_3D)
    xvint = sum(x -> ismissing(x) ? 0.0 : x, x3D .* δz_3D, dims=3) ./ grd.wet3D[:,:,1]
    @series begin
        seriestype := :contourf
        xlabel := "Longitude"
        ylabel := "Latitude"
        colorbar_title := intunit
        lon, lat, view(xvint, :, :, 1)
    end
end

"""
    VerticalAverage(x, grd)

Plots the vertical average of tracer `x`.
"""
@userplot VerticalAverage
@recipe function f(p::VerticalAverage)
    x, grd = p.args
    xunit = string(unit(x[1]))
    x3D = rearrange_into_3Darray(ustrip.(x), grd)
    v = ustrip.(vector_of_volumes(grd))
    v3D = rearrange_into_3Darray(v, grd)
    lon, lat = grd.lon .|> ustrip, grd.lat .|> ustrip
    xmean = sum(x -> ismissing(x) ? 0.0 : x, x3D .* v3D, dims=3) ./
                sum(x -> ismissing(x) ? 0.0 : x, v3D, dims=3)
    @series begin
        seriestype := :contourf
        lon, lat, view(xmean, :, :, 1)
    end
end


"""
    ZonalSlice(x, grd, lon)

Plots a zonal slice of tracer `x` at longitude `lon`.
"""
@userplot ZonalSlice
@recipe function f(p::ZonalSlice)
    x, grd, lon = p.args
    lon = mod(ustrip(lon), 360)
    xunit = string(unit(x[1]))
    x3D = rearrange_into_3Darray(ustrip.(x), grd)
    depth, lat = grd.depth .|> ustrip, grd.lat .|> ustrip
    ix = findfirst(ustrip(grd.lon) .≥ mod(ustrip(lon), 360))
    @series begin
        seriestype := :contourf
        yflip := true
        yticks := Int.(round.(depth))
        ylims := (0, maximum(depth))
        xlabel := "Latitude"
        ylabel := "Depth (m)"
        colorbar_title := xunit
        lat, depth, permutedims(view(x3D,:,ix,:), [2,1])
    end
end

"""
    ZonalAverage(x, grd)

Plots a zonal average of tracer `x`.
"""
@userplot ZonalAverage
@recipe function f(p::ZonalAverage)
    x, grd = p.args
    xunit = string(unit(x[1]))
    x3D = rearrange_into_3Darray(ustrip.(x), grd)
    v = ustrip.(vector_of_volumes(grd))
    v3D = rearrange_into_3Darray(v, grd)
    depth, lat = grd.depth .|> ustrip, grd.lat .|> ustrip
    xmean = sum(x -> ismissing(x) ? 0.0 : x, x3D .* v3D, dims=2) ./
                sum(x -> ismissing(x) ? 0.0 : x, v3D, dims=2)
    @series begin
        seriestype := :contourf
        yflip := true
        yticks := Int.(round.(depth))
        ylims := (0, maximum(depth))
        xlabel := "Latitude"
        ylabel := "Depth (m)"
        colorbar_title := xunit
        lat, depth, permutedims(view(xmean, :, 1, :), [2,1])
    end
end

@userplot ZonalAverage2
@recipe function f(p::ZonalAverage2)
    x, grd = p.args
    x3D = rearrange_into_3Darray(x, grd)
    v = ustrip.(vector_of_volumes(grd))
    v3D = rearrange_into_3Darray(v, grd)
    depth, lat = grd.depth .|> ustrip, grd.lat .|> ustrip
    xmean = sum(x -> ismissing(x) ? 0.0 : x, x3D .* v3D, dims=2) ./
                sum(x -> ismissing(x) ? 0.0 : x, v3D, dims=2)
    xmean = dropdims(xmean, dims=2)

    x, grd, ct = p.args
    x3D = rearrange_into_3Darray(x, grd)
    depths = ustrip.(grd.depth)
    ndepths = length(depths)

    itp = interpolate(xmean, (BSpline(Linear()), BSpline(Linear()))) # interpolate linearly between the data points
    lats = range(ustrip(grd.lat[1]), ustrip(grd.lat[end]), length=length(grd.lat))
    ndepths = length(depth)
    stp = scale(itp, lats, 1:ndepths) # re-scale to the actual domain
    etp = extrapolate(stp, (Line(), Line())) # periodic longitude
    @series begin
        seriestype := :contourf
        yflip := true
        yticks := Int.(round.(depth))
        ylims := (0, maximum(depth))
        y, depth, [etp(lat, i) for i in 1:nz, lat in y]
    end
end

@userplot DepthProfile
@recipe function f(p::DepthProfile)
    x, grd, lat, lon = p.args
    x3D = rearrange_into_3Darray(x, grd)
    depth = ustrip.(grd.depth)
    ix = findfirst(ustrip.(grd.lon) .≥ mod(ustrip(lon), 360))
    iy = findfirst(ustrip.(grd.lat) .≥ ustrip(lat))
    @series begin
        yflip := true
        yticks := Int.(round.(depth))
        ylims := (0, maximum(depth))
        view(x3D, iy, ix, :), depth
    end
end

"""
    InterpolatedDepthProfile(x, grd, lat, lon)

Plots the profile of tracer `x` interpolated at `(lat,lon)` coordinates.
"""
@userplot InterpolatedDepthProfile
@recipe function f(p::InterpolatedDepthProfile)
    x, grd, lat, lon = p.args
    x3D = rearrange_into_3Darray(x, grd)
    depths = ustrip.(grd.depth)
    knots = (ustrip.(grd.lat), ustrip.(grd.lon), 1:length(depths))
    itp = interpolate(knots, x3D, (Gridded(Linear()), Gridded(Linear()), NoInterp()))
    @series begin
        yflip := true
        yticks := Int.(round.(depths))
        ylims := (0, maximum(depths))
        itp(ustrip(lat), ustrip(lon), 1:length(depths)), depths
    end
end


@userplot TransectContourfByDistance
@recipe function f(p::TransectContourfByDistance)
    x, grd, ct = p.args
    x3D = rearrange_into_3Darray(x, grd)
    depths = ustrip.(grd.depth)
    ndepths = length(depths)

    itp = interpolate(periodic_longitude(x3D), (BSpline(Linear()), BSpline(Linear()), NoInterp())) # interpolate linearly between the data points
    lats = range(ustrip(grd.lat[1]), ustrip(grd.lat[end]), length=length(grd.lat))
    lons = range(ustrip(grd.lon[1]), ustrip(grd.lon[1])+360, length=length(grd.lon)+1)
    stp = scale(itp, lats, lons, 1:ndepths) # re-scale to the actual domain
    etp = extrapolate(stp, (Line(), Periodic(), Line())) # periodic longitude

    n = length(ct)
    distances = cumsum(Distances.colwise(Haversine(6371.0),
                [view(ct.lon, :)'; view(ct.lat, :)'],
                [view(ct.lon, [1;1:n-1])'; view(ct.lat, [1;1:n-1])']))
    idx = unique(i -> distances[i], 1:length(ct)) # Remove stations at same (lat,lon)

    @series begin
        seriestype := :contourf
        yflip := true
        yticks := Int.(round.(depths))
        ylims := (0, maximum(depths))
        xlabel := "$(ct.name) distance (km)"
        ylabel := "Depth (m)"
        distances[idx], depths, [etp(lat, lon, i) for i in 1:ndepths, (lat, lon) in zip(ct.lat[idx], ct.lon[idx])]
    end
end

periodic_longitude(x3D::Array{T,3}) where T = view(x3D,:,[1:size(x3D,2); 1],:)
periodic_longitude(x2D::Array{T,2}) where T = view(x2D,:,[1:size(x3D,2); 1])


"""
    ZonalTransect(x, grd, ct)

Plots a zonal transect of tracer `x` along cruise track `ct`.
"""
@userplot ZonalTransect
@recipe function f(p::ZonalTransect)
    x, grd, ct = p.args
    x3D = rearrange_into_3Darray(x, grd)
    depths = ustrip.(grd.depth)
    ndepths = length(depths)

    itp = interpolate(periodic_longitude(x3D), (BSpline(Linear()), BSpline(Linear()), NoInterp())) # interpolate linearly between the data points

    lats = range(ustrip(grd.lat[1]), ustrip(grd.lat[end]), length=length(grd.lat))
    lons = range(ustrip(grd.lon[1]), ustrip(grd.lon[1])+360, length=length(grd.lon)+1)
    stp = scale(itp, lats, lons, 1:ndepths) # re-scale to the actual domain
    etp = extrapolate(stp, (Line(), Periodic(), Line())) # periodic longitude

    isort = sortperm(ct.lat)
    idx = isort[unique(i -> ct.lat[isort][i], 1:length(ct))] # Remove stations at same (lat,lon)
    @series begin
        seriestype := :contourf
        yflip := true
        yticks := Int.(round.(depths))
        ylims := (0, maximum(depths))
        title := "$(ct.name)"
        ylabel := "Depth (m)"
        xlabel := "Latitude (°)"
        ct.lat[idx], depths, [etp(lat, lon, i) for i in 1:ndepths, (lat, lon) in zip(ct.lat[idx], ct.lon[idx])]
    end
end

"""
    MeridionalTransect(x, grd, ct)

Plots a meridional transect of tracer `x` along cruise track `ct`.
"""
@userplot MeridionalTransect
@recipe function f(p::MeridionalTransect)
    x, grd, ct = p.args
    x3D = rearrange_into_3Darray(x, grd)
    depths = ustrip.(grd.depth)
    ndepths = length(depths)

    itp = interpolate(periodic_longitude(x3D), (BSpline(Linear()), BSpline(Linear()), NoInterp())) # interpolate linearly between the data points

    lats = range(ustrip(grd.lat[1]), ustrip(grd.lat[end]), length=length(grd.lat))
    lons = range(ustrip(grd.lon[1]), ustrip(grd.lon[1])+360, length=length(grd.lon)+1)
    stp = scale(itp, lats, lons, 1:ndepths) # re-scale to the actual domain
    etp = extrapolate(stp, (Line(), Periodic(), Line())) # periodic longitude

    shiftedlon = 2sum(90 .≤ ct.lon .≤ 270) > length(ct) ? ct.lon : mod.(ct.lon .+ 180, 360) .- 180

    isort = sortperm(shiftedlon)
    idx = isort[unique(i -> shiftedlon[isort][i], 1:length(ct))] # Remove stations at same (lat,lon)
    @series begin
        seriestype := :contourf
        yflip := true
        yticks := Int.(round.(depths))
        ylims := (0, maximum(depths))
        title := "$(ct.name)"
        ylabel := "Depth (m)"
        xlabel := "Longitude (°)"
        shiftedlon[idx], depths, [etp(lat, lon, i) for i in 1:ndepths, (lat, lon) in zip(ct.lat[idx], shiftedlon[idx])]
    end
end

#=




=#