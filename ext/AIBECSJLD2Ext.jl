module AIBECSJLD2Ext

using AIBECS
using AIBECS: OCIM0, OCIM1, OCIM2, OCCA
using DataDeps
using JLD2
using Unitful

function OCIM0.load()
    OCIM0.register_OCIM0()
    jld2_file = @datadep_str string("AIBECS-OCIM0.1/", "OCIM0.1.jld2")
    @info """You are about to use the OCIM0.1 model.
    If you use it for research, please cite:

    - $(OCIM0.CITATION)

    You can find the corresponding BibTeX entries in the CITATION.bib file
    at the root of the AIBECS.jl package repository.
    (Look for the "DeVries_Primeau_2011" and "Primeau_etal_2013" keys.)
    """
    return jldopen(jld2_file) do file
        file["grid"], ustrip.(file["T"])
    end
end

function OCIM1.load(; version = OCIM1.VERSIONS[1])
    OCIM1.register_OCIM1(; version)
    jld2_file = @datadep_str string("AIBECS-OCIM1_$version/", "OCIM1_$version.jld2")
    @info """You are about to use the OCIM1_$version model.
    If you use it for research, please cite:

    $(OCIM1.CITATION)

    You can find the corresponding BibTeX entries in the CITATION.bib file
    at the root of the AIBECS.jl package repository.
    (Look for the "DeVries_Primeau_2011" and "DeVries_2014" keys.)
    """
    return jldopen(jld2_file) do file
        file["grid"], ustrip.(file["T"])
    end
end

function OCIM2.load(; version = OCIM2.VERSIONS[1], HeFluxes = false)
    OCIM2.register_OCIM2(; version)
    jld2_file = @datadep_str string("AIBECS-OCIM2_$version/", "OCIM2_$version.jld2")
    @info """You are about to use the OCIM2_$version model.
    If you use it for research, please cite:

    $(OCIM2.CITATION)

    You can find the corresponding BibTeX entries in the CITATION.bib file
    at the root of the AIBECS.jl package repository.
    (Look for the "DeVries_Holzer_2019" key.)
    """
    return if HeFluxes
        jldopen(jld2_file) do file
            file["grid"], ustrip.(file["T"]), file["He3Flux"], file["He4Flux"]
        end
    else
        jldopen(jld2_file) do file
            file["grid"], ustrip.(file["T"])
        end
    end
end

function OCCA.load()
    OCCA.register_OCCA()
    jld2_file = @datadep_str string("AIBECS-OCCA/", "OCCA.jld2")
    @info """You are about to use the OCCA model.
    If you use it for research, please cite:

    - $(OCCA.CITATION)

    You can find the corresponding BibTeX entries in the CITATION.bib file
    at the root of the AIBECS.jl package repository.
    (Look for the "Forget_2010" key.)
    """
    return jldopen(jld2_file) do file
        file["grid"], ustrip.(file["T"])
    end
end

end # module
