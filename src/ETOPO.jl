module ETOPO

using OceanGrids
using Dates
using Unitful
using DataDeps              # For storage location of data
using NCDatasets
using Distances

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Base.download(remotepath, localpath)
    return localpath
end

# Create registry entry for ETOPO
function register_ETOPO()
    register(
        DataDep(
            "ETOPO",
            """
            Reference:
            $(citation()) [accessed $(Dates.today())]
            """,
            "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gdal.grd.gz",
            sha2_256,
            fetch_method = fallback_download,
            post_fetch_method = unpack
        )
    )
    return nothing
end

"""
    Z, lats, lons = load()

Returns the fine resolution topography from ETOPO.
"""
function load()
    register_ETOPO()
    nc_file = @datadep_str string("ETOPO/ETOPO1_Bed_g_gdal.grd")
    ds = Dataset(nc_file)
    Z, lats, lons, dims =
    Dataset(nc_file,"r") do ds
        dims = Tuple(ds["dimension"][:])
        ds["z"][:], range(ds["y_range"]..., length=dims[2]),
        range(ds["x_range"]..., length=dims[1]), dims
    end # ds is closed
    @info """You are about to use the ETOPO data set.
          If you use it for research, please cite:

          - $(citation())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "Amante_Eakins_2009" key.)
          """
    return reverse(reshape(Z, dims), dims=2), lats, lons # TODO
end

citation() = "Amante, C. and B.W. Eakins, 2009. ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis. NOAA Technical Memorandum NESDIS NGDC-24. National Geophysical Data Center, NOAA. doi:10.7289/V5C8276M [access date]."

"""
    bintopo(Z, lat, lon, grd)

Bins the depths `Z` into at locations `(lat,lon)` onto `grd`.
To be used with the ETOPO dataset to dertermine the amount of subgrid topography in each box area
"""
function bintopo(Z::AbstractArray{T,2}, lats, lons, grd) where T
    bin3D = zeros(size(grd))
    # Use the top edge to invert the functions lon[i], lat[j], and depth[k]
    # I.e., to find the index in grd for a given lon, lat, depth
    edepth = ustrip.(grd.depth .+ 0.5grd.δdepth)
    elon = ustrip.(grd.lon .+ 0.5grd.δlon)
    elat = ustrip.(grd.lat .+ 0.5grd.δlat)
    kmax = dropdims(count(iswet(grd), dims=3), dims=3)
    for (j,lat) in enumerate(lats)
        jgrd = searchsortedfirst(elat, lat)
        for (i,lon) in enumerate(lons)
            igrd = searchsortedfirst(elon, mod(lon, 360))
            kmax[jgrd,igrd] == 0 && continue
            kgrd = searchsortedfirst(view(edepth, 1:kmax[jgrd, igrd]-1), -Z[i,j])
            bin3D[jgrd, igrd, kgrd] += 1
        end
    end
    return bin3D
end

"""
    fractiontopo(grd)

Returns the fraction of subgrid sediments in each wet box of grd.

Note that (i) points located above sea level are counted in the top layer,
and (ii) points located in a dry box or below are counted in the deepest
wet box.
"""
function fractiontopo(grd)
    @warn "Binning ETOPO to grd. This may take a few seconds"
    Z, lats, lons = load()
    return fractiontopo(bintopo(Z, lats, lons, grd), grd)
end
fractiontopo(bin3D, grd) = vectorize(bin3D ./ repeat(sum(bin3D, dims=3), outer=(1, 1, size(grd)[3])), grd)
export fractiontopo

function roughnesstopo(grd)
    @warn "Estimating roughness from ETOPO. This may take a few seconds"
    Z, lats, lons = load()
    R_earth = ustrip(u"m", 6371.0u"km")
    δlon = step(lons)
    δlat = step(lons)
    midlats = 0.5(lats[1:end-1] + lats[2:end])
    midlons = 0.5(lons[1:end-1] + lons[2:end])
    δxs = [Haversine(R_earth)((0.0, lat), (δlon, lat)) for lat in midlats]
    δy = Haversine(R_earth)((lons[1], 0.0), (lons[1], δlat))
    A_lon = abs.(diff(Z, dims=1)) .* δy
    Z_lon = 0.5(Z[1:end-1,:] + Z[2:end,:])
    A_lat = abs.(diff(Z, dims=2)) .* repeat(reshape(δxs, (1, length(δxs))), outer=(size(Z,1), 1))
    Z_lat = 0.5(Z[:,1:end-1] + Z[:,2:end])
    return sumroughness(A_lon, Z_lon, lats, midlons, grd) + sumroughness(A_lat, Z_lat, midlats, lons, grd)
end

function sumroughness(A::AbstractArray{T,2}, Z, lats, lons, grd) where T
    roughness = zeros(size(grd))
    edepth = ustrip.(grd.depth .+ 0.5grd.δdepth)
    elon = ustrip.(grd.lon .+ 0.5grd.δlon)
    elat = ustrip.(grd.lat .+ 0.5grd.δlat)
    kmax = dropdims(count(iswet(grd), dims=3), dims=3)
    for (j,lat) in enumerate(lats)
        jgrd = searchsortedfirst(elat, lat)
        for (i,lon) in enumerate(lons)
            igrd = searchsortedfirst(elon, mod(lon, 360))
            kmax[jgrd,igrd] == 0 && continue
            kgrd = searchsortedfirst(view(edepth, 1:kmax[jgrd, igrd]-1), -Z[i,j])
            roughness[jgrd, igrd, kgrd] += A[i,j]
        end
    end
    return roughness[iswet(grd)]
end

end # module

export ETOPO