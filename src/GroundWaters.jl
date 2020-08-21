module GroundWaters

using OceanGrids
using Unitful
using Unitful: m, yr
using DataDeps              # For storage location of data
using Shapefile
using DataFrames
import OceanGrids: regrid

struct GroundWaterSource{T}
    lon::Float64
    lat::Float64
    VFR::T # volumetric flow rate
end

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Base.download(remotepath, localpath)
    return localpath
end

# Create registry entry for ETOPO
function register_groundwater_discharge()
    register(
        DataDep(
            "groundwater_discharge",
            """
            Reference:
            $(citation())
            """,
            "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-15064-8/MediaObjects/41467_2020_15064_MOESM8_ESM.zip",
            "aed5141cf00cfd6c2f7f70d4a18582b4613de571737661287ca09785c16f7e05",
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
    register_groundwater_discharge()
    shp_file = @datadep_str string("groundwater_discharge/Supplementary_data_3/coastal_gw_discharge.shp")
    df = DataFrame(Shapefile.Table(shp_file))
    gws = GroundWaterSource.(df[:, :ws_cent_x], df[:, :ws_cent_y], df[:, :cgd3_best] * m^3/yr)
    gws = [x for x in gws if x.VFR ≥ 0m^3/yr]
    @info """You are about to use the groundwater discharge data set.
          If you use it for research, please cite:

          $(citation())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "Luijendijk_etal_2019" and "Luijendijk_etal_2020" keys.)
          """
    return gws
end

citation() = """
- Luijendijk, E., Gleeson, T. & Moosdorf, N. Fresh groundwater discharge insignificant for the world’s oceans but important for coastal ecosystems. Nat Commun 11, 1260 (2020). doi:10.1038/s41467-020-15064-8
- Luijendijk, Elco; Gleeson, Tom; Moosdorf, Nils (2019): Geospatial data and model results for a global model study of coastal groundwater discharge. PANGAEA, doi:10.1594/PANGAEA.907641
"""

"""
    regrid(gws, grd)

Regrids and bins the groundwater discharge values into `grd` surface boxes.
"""
function regrid(gws::Vector{GroundWaterSource{T}}, grd) where T <: Quantity
    lats = [x.lat for x in gws]
    lons = [x.lon for x in gws]
    depths = zeros(length(gws))
    vs = [x.VFR for x in gws]
    return regrid(vs, lats, lons, depths, grd)
end
export regrid

end # module

export GroundWaters