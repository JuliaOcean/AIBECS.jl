module AIBECSJLD2Ext

using AIBECS
using AIBECS: OCIM0, OCIM1, OCIM2, OCCA
using DataDeps
using JLD2
using LinearAlgebra: Diagonal
using Unitful

function OCIM0.load()
    OCIM0.register_OCIM0()
    OCIM0.invalidate_stale_cache()
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
    OCIM1.invalidate_stale_cache(; version)
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
    OCIM2.invalidate_stale_cache(; version)
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

function OCCA.load(; conservemass = true)
    OCCA.register_OCCA()
    OCCA.invalidate_stale_cache()
    jld2_file = @datadep_str string("AIBECS-OCCA/", "OCCA.jld2")
    @info """You are about to use the OCCA model.
    If you use it for research, please cite:

    - $(OCCA.CITATION)

    You can find the corresponding BibTeX entries in the CITATION.bib file
    at the root of the AIBECS.jl package repository.
    (Look for the "Forget_2010" key.)
    """
    grd, T = jldopen(jld2_file) do file
        file["grid"], ustrip.(file["T"])
    end
    if conservemass
        @warn "Enforcing mass conservation: modifying OCCA transport matrix as T -= Diagonal((Tᵀ v) ./ v). Pass `conservemass = false` to skip."
        v = volumevec(grd)
        T = T - Diagonal((T' * v) ./ v)
    else
        @warn "You are about to use the raw OCCA transport matrix, which does not conserve mass. Pass `conservemass = true` (which is the default) to enforce mass conservation."
    end
    return grd, T
end

end # module
