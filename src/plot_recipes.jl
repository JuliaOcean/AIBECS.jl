


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
        seriestype := :heatmap
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
    xguide --> "Longitude"
    yguide --> "Latitude"
    @series begin
        seriestype := :heatmap
        grd.lon, grd.lat, view(x3D, :, :, 1)
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
    x2D = nansum(x3D .* δz_3D, dims=3) ./ grd.wet3D[:,:,iz1] # 0/0 = NaN, 0/1=0, x/0=Inf, and x/1=x
    @series begin
        seriestype := :heatmap
        xguide --> "Longitude"
        yguide --> "Latitude"
        grd.lon, grd.lat, view(x2D, :, :, 1)
    end
end
finddepthindex(depths::Tuple, grd) = ((finddepthindex(d, grd) for d in depths)...,)
convertdepth(depths::Tuple) = ((convertdepth(d) for d in depths)...,)

"""
    VerticalSum(x, grd [; depthlim])

Plots the vertical integral of tracer `x`.
"""
@userplot VerticalSpecialSum
@recipe function f(p::VerticalSpecialSum; depthlim=extrema(p.args[2].depth))
    x, grd = p.args
    depthlim = convertdepth(depthlim)
    iz1, iz2 = finddepthindex(depthlim, grd)
    x3D = rearrange_into_3Darray(x, grd)[:,:,iz1:iz2]
    v = vector_of_volumes(grd)
    v3D = rearrange_into_3Darray(v, grd)[:,:,iz1:iz2]
    δz_3D = view(grd.δz_3D, :, :, iz1:iz2)
    x2D = nansum(x3D .* δz_3D, dims=3) ./ nansum(δz_3D .* grd.wet3D[:,:,iz1:iz2], dims=3)
    @series begin
        seriestype := :heatmap
        xguide --> "Longitude"
        yguide --> "Latitude"
        grd.lon, grd.lat, view(x2D, :, :, 1)
    end
end

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
        seriestype := :heatmap
        xguide --> "Longitude"
        yguide --> "Latitude"
        grd.lon, grd.lat, view(x2D, :, :, 1)
    end
end

# sum skipping nans
nansum(x; kwargs...) = sum(x .* (!isnan).(x); kwargs...)

"""
    MeridionalSlice(x, grd; lon)

Plots a zonal slice of tracer `x` at longitude `lon`.
"""
@userplot MeridionalSlice
@recipe function f(p::MeridionalSlice; lon=nothing)
    isnothing(lon) && error("you must specify the longitude, e.g., `lon=100`")
    x, grd = p.args
    lon = convertlon(lon)
    ix = findlonindex(lon, grd)
    x3D = rearrange_into_3Darray(x, grd)
    @series begin
        seriestype := :heatmap
        yflip := true
        ylims --> (0, 1) .* sum(grd.δdepth)
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
    ZonalSlice(x, grd, lat)

Plots a Meridional slice of tracer `x` at longitude `lat`.
"""
@userplot ZonalSlice
@recipe function f(p::ZonalSlice; lat=nothing)
    isnothing(lat) && error("you must specify the latitude, e.g., `lat=-30`")
    x, grd = p.args
    lat = convertlat(lat)
    iy = findlatindex(lat, grd)
    x3D = rearrange_into_3Darray(x, grd)
    @series begin
        seriestype := :heatmap
        yflip := true
        ylims --> (0, 1) .* sum(grd.δdepth)
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
        seriestype := :heatmap
        yflip := true
        ylims --> (0, 1) .* sum(grd.δdepth)
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
        ylims --> (0, 1) .* sum(grd.δdepth)
        yguide --> "Depth"
        itp(lat, lon, 1:ndepth), grd.depth
    end
end

# TODO fix this plot to smarlty reorder contour.
@userplot TransectByDistance
@recipe function f(p::TransectByDistance; ct=nothing, start=:south, show_stations=false)
    isnothing(ct) && error("you must include the cruise track with `ct=...`")
    x, grd = p.args
    x3D = lonextend(rearrange_into_3Darray(x, grd))
    knots = ustrip.((grd.lat, lonextend(grd.lon), grd.depth))
    itp = interpolate(knots, x3D, Gridded(Constant()))

    ct = sort(ct, start=start)
    ctlats = [convertlat(st.lat) for st in ct.stations]
    ctlons = [uconvertlon(st.lon) for st in ct.stations]

    n = length(ct)
    distance = vcat(0.0u"km", cumsum(diff(ct)))
    ikeep = findall(diff(distance) .> 0u"km")
    push!(ikeep, n)
    distance = distance[ikeep]
    ctlats = ctlats[ikeep]
    ctlons = ctlons[ikeep]
    nkeep = length(ikeep)
    itplats = LinearInterpolation(1:nkeep, ctlats)
    itplons = LinearInterpolation(1:nkeep, ctlons)
    itpdist = LinearInterpolation(1:nkeep, distance)
    finedist = itpdist(range(1, nkeep, length=10nkeep))
    finelats = itplats(range(1, nkeep, length=10nkeep))
    finelons = itplons(range(1, nkeep, length=10nkeep))

    yflip := true
    title --> "$(ct.name)"
    yguide --> "Depth"
    xguide --> "Distance"
    if show_stations
        @series begin
            markershape := :vline
            label := ""
            color := :black
            distance, zeros(nkeep)
        end
    end
    @series begin
        seriestype := :heatmap
        yflip := true
        ylims --> (0, 1) .* sum(grd.δdepth)
        finedist, grd.depth, [itp(ustrip.((lat, lon, d))...) for d in grd.depth, (lat,lon) in zip(finelats, finelons)]
    end
    if haskey(plotattributes, :levels)
        @series begin
            seriestype := :contour
            color := :black
            contour_labels := true
            levels := plotattributes[:levels]
            distance, grd.depth, [itp(ustrip.((lat, lon, d))...) for d in grd.depth, (lat,lon) in zip(ctlats, ctlons)]
        end    
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
    t = p.args[1]
    if hasproperty(t, :stations)
        lons = [st.lon for st in t.stations]
        lats = [st.lat for st in t.stations]
    elseif hasproperty(t, :profiles)
        lons = [pro.station.lon for pro in t.profiles]
        lats = [pro.station.lat for pro in t.profiles]
    end
    xguide --> "Longitude"
    yguide --> "Latitude"
    @series begin
        label --> (hasproperty(t, :name) ? t.name : (hasproperty(t, :cruise) ? t.cruise : nothing))
        linewidth --> 3
        linecolor --> :black
        uconvertlon.(lons), uconvertlat.(lats)
    end
    @series begin
        label --> nothing
        markershape --> [:circle, :circle]
        markercolor --> [:green, :red]
        linewidth --> 0
        uconvertlon.(lons)[[1,end]], uconvertlat.(lats)[[1,end]]
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
    x3D = lonextend(rearrange_into_3Darray(x, grd))
    knots = ustrip.((grd.lat, lonextend(grd.lon), grd.depth))
    itp = interpolate(knots, x3D, Gridded(Constant()))
    ctlats = [convertlat(st.lat) for st in ct.stations]
    ctlons = [uconvertlon(st.lon) for st in ct.stations]
    isort = sortperm(ctlats)
    idx = isort[unique(i -> ctlats[isort][i], 1:length(ct))] # Remove stations at same (lat,lon)
    ctlats2 = ctlats[idx]
    ctlons2 = ctlons[idx]
    itplons = LinearInterpolation(ctlats2, ctlons2)
    finelats = range(ctlats2[1], ctlats2[end], length=100length(ctlats2))
    @series begin
        seriestype := :heatmap
        yflip := true
        ylims --> (0, 1) .* sum(grd.δdepth)
        title --> "$(ct.name)"
        yguide --> "Depth"
        xguide --> "Latitude"
        finelats, grd.depth, [itp(ustrip.((lat, itplons(lat), d))...) for d in grd.depth, lat in finelats]
    end
end
uconvertlon(lon::T where {T<:Real}) = lon * u"°"
uconvertlon(lon::T where {T<:Quantity}) = unit(lon) == Unitful.° ? lon : error("Not a valid lon")
uconvertlat(lat::T where {T<:Real}) = lat * u"°"
uconvertlat(lat::T where {T<:Quantity}) = unit(lat) == Unitful.° ? lat : error("Not a valid lat")

function lonextend(x3D::Array{T,N} where {T,N})
    x3Dext = Array{eltype(x3D),3}(undef, (size(x3D) .+ (0,2,0))...)
    x3Dext[:,2:end-1,:] .= x3D
    x3Dext[:,1,:] .= x3D[:,end,:]
    x3Dext[:,end,:] .= x3D[:,1,:]
    return x3Dext
end
function lonextend(lon::Vector{T} where {T})
    lonext = Vector{eltype(lon)}(undef, length(lon) + 2)
    lonext[2:end-1] .= lon
    lonext[1] = lon[end] - 360u"°"
    lonext[end] = lon[1] + 360u"°"
    return lonext
end

"""
    MeridionalScatterTransect(t)

Plots a scatter of the discrete obs of `t` in (lat,depth) space.
"""
@userplot MeridionalScatterTransect
@recipe function f(p::MeridionalScatterTransect)
    transect = p.args[1]
    depths = reduce(vcat, pro.depths for pro in transect.profiles)
    values = reduce(vcat, pro.values for pro in transect.profiles)
    lats = reduce(vcat, pro.station.lat * ones(length(pro)) for pro in transect.profiles)
    @series begin
        seriestype := :scatter
        yflip := true
        marker_z --> values
        markershape --> :circle
        xlims --> extrema(lats)
        clims --> (0, 1) .* maximum(transect)
        label --> "$(transect.tracer) along $(transect.cruise)"
        yguide --> "Depth"
        xguide --> "Latitude"
        convertlat.(lats), convertdepth.(depths)
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
    values = reduce(vcat, pro.values for pro in transect.profiles)
    lons = reduce(vcat, pro.station.lon * ones(length(pro)) for pro in transect.profiles)
    @series begin
        seriestype := :scatter
        yflip := true
        marker_z --> values
        markershape --> :circle
        xlims --> extrema(lons)
        clims --> (0, 1) .* maximum(transect)
        label --> "$(transect.tracer) along $(transect.cruise)"
        yguide --> "Depth"
        xguide --> "Longitude"
        uconvertlon.(lons), convertdepth.(depths)
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
    x3D = lonextend(rearrange_into_3Darray(x, grd))
    knots = ustrip.((grd.lat, lonextend(grd.lon), grd.depth))
    itp = interpolate(knots, x3D, Gridded(Constant()))
    ctlats = [convertlat(st.lat) for st in ct.stations]
    ctlons = [uconvertlon(st.lon) for st in ct.stations]
    isort = sortperm(ctlons)
    idx = isort[unique(i -> ctlons[isort][i], 1:length(ct))] # Remove stations at same (lat,lon)
    ctlats2 = ctlats[idx]
    ctlons2 = ctlons[idx]
    itplats = LinearInterpolation(ctlons2, ctlats2)
    finelons = range(ctlons2[1], ctlons2[end], length=100length(ctlons2))
    finelats = itplats.(finelons)
    @series begin
        seriestype := :heatmap
        yflip := true
        ylims --> (0, 1) .* sum(grd.δdepth)
        title --> "$(ct.name)"
        yguide --> "Depth"
        xguide --> "Longitude"
        finelons, grd.depth, [itp(ustrip.((itplats(lon), mod(lon, 360u"°"), d))...) for d in grd.depth, lon in finelons]
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