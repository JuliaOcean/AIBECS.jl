module AIBECSETOPOExt

using AIBECS
using AIBECS.ETOPO: ETOPO, register_ETOPO, citation, DATA_FILE_NAME, bintopo, sumroughness
using DataDeps
using Distances: Haversine
using NCDatasets
using Unitful
using OceanGrids: iswet

function ETOPO.load(data = :bedrock)
    register_ETOPO(data)
    nc_file = @datadep_str joinpath(string("ETOPO_", data), DATA_FILE_NAME[data])
    Z, lats, lons, dims = Dataset(nc_file, "r") do ds
        dims = Tuple(ds["dimension"][:])
        ds["z"][:], range(ds["y_range"]..., length = dims[2]),
            range(ds["x_range"]..., length = dims[1]), dims
    end
    @info """You are about to use the ETOPO data set.
    If you use it for research, please cite:

    - $(citation())

    You can find the corresponding BibTeX entries in the CITATION.bib file
    at the root of the AIBECS.jl package repository.
    (Look for the "Amante_Eakins_2009" key.)
    """
    return reverse(reshape(Z, dims), dims = 2), lats, lons
end

function ETOPO.fractiontopo(grd)
    @warn "Binning ETOPO to grd. This may take a few seconds"
    Z, lats, lons = ETOPO.load()
    return ETOPO.fractiontopo(bintopo(Z, lats, lons, grd), grd)
end

function ETOPO.roughnesstopo(grd)
    @warn "Estimating roughness from ETOPO. This may take a few seconds"
    Z, lats, lons = ETOPO.load()
    R_earth = ustrip(u"m", 6371.0u"km")
    δlon = step(lons)
    δlat = step(lons)
    midlats = 0.5(lats[1:(end - 1)] + lats[2:end])
    midlons = 0.5(lons[1:(end - 1)] + lons[2:end])
    δxs = [Haversine(R_earth)((0.0, lat), (δlon, lat)) for lat in midlats]
    δy = Haversine(R_earth)((lons[1], 0.0), (lons[1], δlat))
    A_lon = abs.(diff(Z, dims = 1)) .* δy
    Z_lon = 0.5(Z[1:(end - 1), :] + Z[2:end, :])
    A_lat = abs.(diff(Z, dims = 2)) .* repeat(reshape(δxs, (1, length(δxs))), outer = (size(Z, 1), 1))
    Z_lat = 0.5(Z[:, 1:(end - 1)] + Z[:, 2:end])
    return sumroughness(A_lon, Z_lon, lats, midlons, grd) +
        sumroughness(A_lat, Z_lat, midlats, lons, grd)
end

end # module
