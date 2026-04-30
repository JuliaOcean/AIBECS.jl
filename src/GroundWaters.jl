"""
    GroundWaters

Coastal fresh-groundwater discharge dataset of Luijendijk et al. (2020),
exposed as `GroundWaterSource` entries that can be regridded onto an
`OceanGrid`. Loading the dataset requires the `Shapefile` and `DataFrames`
packages so that the `AIBECSShapefileExt` extension is activated.
"""
module GroundWaters

using OceanGrids
using Unitful
using Unitful: m, yr
using DataDeps              # For storage location of data
using Downloads
import OceanGrids: regrid

"""
    GroundWaterSource(lon, lat, VFR)

A point-source groundwater discharge: longitude, latitude, and volumetric flow
rate `VFR` (Unitful quantity).
"""
struct GroundWaterSource{T}
    lon::Float64
    lat::Float64
    VFR::T # volumetric flow rate
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
            "aed5141cf00cfd6c2f7f70d4a18582b4613de571737661287ca09785c16f7e05";
            post_fetch_method = unpack
        )
    )
    return nothing
end

"""
    gws = load()

Returns the coastal groundwater discharge dataset.

Requires `using Shapefile, DataFrames` so that the `AIBECSShapefileExt`
extension is activated.
"""
function load(args...; kwargs...)
    error("AIBECS.GroundWaters.load requires `using Shapefile, DataFrames`. " *
          "Add them to your environment, then retry.")
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