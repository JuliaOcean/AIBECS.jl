"""
    ETOPO

Access to the 1-arc-minute ETOPO1 global relief model (Amante and Eakins, 2009),
either the bedrock (`:bedrock`) or ice-surface (`:ice`) variant. Provides
[`ETOPO.load`](@ref) plus the bin-on-grid helpers [`fractiontopo`](@ref) and
[`roughnesstopo`](@ref) for mapping subgrid topography onto an `OceanGrid`.

Loading the dataset requires `Distances` and `NCDatasets` so that the
`AIBECSETOPOExt` extension is activated.
"""
module ETOPO

using OceanGrids
using Dates
using Unitful
using DataDeps              # For storage location of data
using Downloads

const DATA_FILE_NAME = Dict(
    :bedrock => "ETOPO1_Bed_g_gdal.grd",
    :ice => "ETOPO1_Ice_g_gdal.grd"
)

const URL = Dict(
    :bedrock => "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gdal.grd.gz",
    :ice => "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gdal.grd.gz"
)

# Create registry entry for ETOPO
function register_ETOPO(data=:bedrock)
    register(
        DataDep(
            string("ETOPO_", data),
            """
            Reference:
            $(citation()) [accessed $(Dates.today())]
            """,
            URL[data],
            sha2_256;
            post_fetch_method = unpack
        )
    )
    return nothing
end


"""
    Z, lats, lons = ETOPO.load()
    Z, lats, lons = ETOPO.load(:bedrock)
    Z, lats, lons = ETOPO.load(:ice)

Returns the fine resolution topography from ETOPO.

Requires `using Distances, NCDatasets` so that the `AIBECSETOPOExt`
extension is activated.
"""
function load(args...; kwargs...)
    error("AIBECS.ETOPO.load requires `using Distances, NCDatasets`. " *
          "Add them to your environment, then retry.")
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

Requires `using Distances, NCDatasets` so that the `AIBECSETOPOExt`
extension is activated. Note that (i) points located above sea level
are counted in the top layer, and (ii) points located in a dry box or
below are counted in the deepest wet box.
"""
function fractiontopo end
# Pure-math helper: takes the binned ETOPO array and returns the fraction
# per wet box. No NCDatasets/Distances dependency, so it stays in main.
fractiontopo(bin3D, grd) = vectorize(bin3D ./ repeat(sum(bin3D, dims=3), outer=(1, 1, size(grd)[3])), grd)
export fractiontopo

"""
    roughnesstopo(grd)

Estimates topographic roughness from ETOPO.

Requires `using Distances, NCDatasets` so that the `AIBECSETOPOExt`
extension is activated.
"""
function roughnesstopo end

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