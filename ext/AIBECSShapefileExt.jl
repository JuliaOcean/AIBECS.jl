module AIBECSShapefileExt

using AIBECS
using AIBECS.GroundWaters: GroundWaters, GroundWaterSource, register_groundwater_discharge, citation
using DataDeps
using Shapefile
using DataFrames
using Unitful
using Unitful: m, yr

function GroundWaters.load()
    register_groundwater_discharge()
    shp_file = @datadep_str string("groundwater_discharge/Supplementary_data_3/coastal_gw_discharge.shp")
    df = DataFrame(Shapefile.Table(shp_file))
    gws = GroundWaterSource.(df[:, :ws_cent_x], df[:, :ws_cent_y], df[:, :cgd3_best] * m^3 / yr)
    gws = [x for x in gws if x.VFR ≥ 0m^3 / yr]
    @info """You are about to use the groundwater discharge data set.
    If you use it for research, please cite:

    $(citation())

    You can find the corresponding BibTeX entries in the CITATION.bib file
    at the root of the AIBECS.jl package repository.
    (Look for the "Luijendijk_etal_2019" and "Luijendijk_etal_2020" keys.)
    """
    return gws
end

end # module
