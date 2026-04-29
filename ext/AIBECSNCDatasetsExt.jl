module AIBECSNCDatasetsExt

using AIBECS
using AIBECS.AeolianSources: AeolianSources, register_AeolianSource,
    DATASET_LONGNAMES, Chien_AEROSOLTYPE_NAMES, Kok_REGIONS_NAMES, CITATIONS
using DataDeps
using NCDatasets

function AeolianSources.load_Chien()
    s = "Chien"
    register_AeolianSource(s)
    nc_file = @datadep_str joinpath("AIBECS-$(DATASET_LONGNAMES[s])", "post.aerosols.2x2.seasonal.nc")
    s_A_2D = Dataset(nc_file, "r") do ds
        Dict(
            :lat => ds["lat"][:],
            :lon => ds["lon"][:],
            (t => ds["dep"][:,:,i] for (i,t) in enumerate(Chien_AEROSOLTYPE_NAMES))...
        )
    end
    @info """You are about to use the Chien et al. (2016) data for aeolian deposition.
          If you use it for research, please cite:

          $(CITATIONS[s])

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the $(DATASET_LONGNAMES[s]) key.)
          """
    return s_A_2D
end

function AeolianSources.load_Kok()
    s = "Kok"
    register_AeolianSource(s)
    nc_file = @datadep_str joinpath("AIBECS-$(DATASET_LONGNAMES[s])", "DustCOMM_source_region_wetdep_annual_PM20_abs.nc")
    s_A_2D = Dataset(nc_file, "r") do ds
        Dict(
            :lat => ds["lat"][:],
            :lon => ds["lon"][:],
            (r => ds["Mean"][i,:,:] for (i,r) in enumerate(Kok_REGIONS_NAMES))...
        )
    end
    @info """You are about to use the Kok et al. (2021) data for aeolian deposition.
          If you use it for research, please cite:

          $(CITATIONS[s])

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the $(DATASET_LONGNAMES[s]) key.)
          """
    return s_A_2D
end

end # module
