@userplot HorizontalSlice
@recipe function f(p::HorizontalSlice)
    x, grd, depth = p.args
    x3D = rearrange_into_3Darray(x, grd)
    lon, lat = grd.lon .|> ustrip, grd.lat .|> ustrip
    iz = findfirst(ustrip.(grd.depth) .≥ ustrip(depth))
    isnothing(iz) && (iz = length(grd.depth)) 
    @series begin
        seriestype := :contourf
        lon, lat, view(x3D, :, :, iz)
    end
end


@userplot SurfaceMap
@recipe function f(p::SurfaceMap)
    x, grd = p.args
    x3D = rearrange_into_3Darray(x, grd)
    lon, lat = grd.lon .|> ustrip, grd.lat .|> ustrip
    @series begin
        seriestype := :contourf
        lon, lat, view(x3D, :, :, 1)
    end
end


@userplot VerticalAverage
@recipe function f(p::VerticalAverage)
    x, grd = p.args
    x3D = rearrange_into_3Darray(x, grd)
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


@userplot ZonalSlice
@recipe function f(p::ZonalSlice)
    x, grd, lon = p.args
    x3D = rearrange_into_3Darray(x, grd)
    depth, lat = grd.depth .|> ustrip, grd.lat .|> ustrip
    ix = findfirst(ustrip(grd.lon) .≥ mod(ustrip(lon), 360))
    @series begin
        seriestype := :contourf
        yflip := true
        yticks := Int.(round.(depth))
        ylims := (0, maximum(depth))
        lat, depth, permutedims(view(x3D,:,ix,:), [2,1])
    end
end


@userplot ZonalAverage
@recipe function f(p::ZonalAverage)
    x, grd = p.args
    x3D = rearrange_into_3Darray(x, grd)
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
        lat, depth, permutedims(view(xmean, :, 1, :), [2,1])
    end
end


@userplot DepthProfile
@recipe function f(p::DepthProfile)
    x, grd, lat, lon = p.args
    x3D = rearrange_into_3Darray(x, grd)
    depth = ustrip.(grd.depth)
    ix = findfirst(ustrip(grd.lon) .≥ mod(ustrip(lon), 360))
    iy = findfirst(ustrip(grd.lat) .≥ ustrip(lat))
    @series begin
        yflip := true
        yticks := Int.(round.(depth))
        ylims := (0, maximum(depth))
        view(x3D, iy, ix, :), depth
    end
end

