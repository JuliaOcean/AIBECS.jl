"""
This `OCIM2` module is used to load the OCIM2 matrices and grid for use in AIBECS.

!!! tip
    To load the default OCIM2 matrix and grid, do
    ```
    julia> grd, T = OCIM2.load()
    ```
    But you can also load the other matrices by specifying which version you want, e.g.,
    ```
    julia> grd, T = OCIM2.load(version="OCIM2_KiHIGH_noHe")
    ```
    See *DeVries and Holzer* (2019) for more details


!!! note
    The files, that are downloaded from a public and persistant URL in FigShare,
    were created with the code available at https://github.com/briochemc/OceanCirculations.
"""
module OCIM2

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using Downloads
using JLD2                  # For saving circulation as JLD2 format
using CodecZlib             # for JLD2 compression JLD2
using Unitful               # for units
using Reexport
@reexport using OceanGrids            # To store the grid

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Downloads.download(remotepath, localpath)
    return localpath
end

# OCIM2 URLs
function url(; version="CTL_He")
    url_start = "https://files.figshare.com/"
    version == "CTL_He"             ? "$(url_start)28336284/OCIM2_CTL_He.jld2" :
    version == "CTL_noHe"           ? "$(url_start)28336299/OCIM2_CTL_noHe.jld2" :
    version == "KiHIGH_He"          ? "$(url_start)28336302/OCIM2_KiHIGH_He.jld2" :
    version == "KiHIGH_noHe"        ? "$(url_start)28336311/OCIM2_KiHIGH_noHe.jld2" :
    version == "KiLOW_He"           ? "$(url_start)28336317/OCIM2_KiLOW_He.jld2" :
    version == "KiLOW_noHe"         ? "$(url_start)28336323/OCIM2_KiLOW_noHe.jld2" :
    version == "KvHIGH_He"          ? "$(url_start)28336326/OCIM2_KvHIGH_He.jld2" :
    version == "KvHIGH_KiHIGH_noHe" ? "$(url_start)28336329/OCIM2_KvHIGH_KiHIGH_noHe.jld2" :
    version == "KvHIGH_KiLOW_He"    ? "$(url_start)28336332/OCIM2_KvHIGH_KiLOW_He.jld2" :
    version == "KvHIGH_KiLOW_noHe"  ? "$(url_start)28336341/OCIM2_KvHIGH_KiLOW_noHe.jld2" :
    version == "KvHIGH_noHe"        ? "$(url_start)28336353/OCIM2_KvHIGH_noHe.jld2" :
    OCIM2versionerror(version)
end

OCIM2versionerror(version) = error("""`$version` is not a valid OCIM2 version name.

                                   Valid versions are `CTL_He`, `CTL_noHe`, `KiHIGH_He`, `KiHIGH_noHe`, `KiLOW_He`, `KiLOW_noHe`, `KvHIGH_He`, `KvHIGH_KiHIGH_noHe`, `KvHIGH_KiLOW_He`, `KvHIGH_KiLOW_noHe`, and `KvHIGH_noHe`.

                                   See *DeVries and Holzer* (2019) for more details on OCIM2 configurations.""")

# OCIM2 Hashes
function sha(; version="CTL_He")
    version == "CTL_He"             ? "610d6abd269950766b75f48ceeee84eb6461782f88331931233bc8da9b96f69d" :
    version == "CTL_noHe"           ? "87c486ae59c2e7702fe2812004e49164ee83dccf00cc05868de0127fc09af3da" :
    version == "KiHIGH_He"          ? "1a2c115696b88df0191126102cea3ba201042e756ec7ed73d1100f0eddefa6c6" :
    version == "KiHIGH_noHe"        ? "cd6e27786f78dd05497fd2c4eb2fb6d94d715940c81378ecefa3c0da35e2977b" :
    version == "KiLOW_He"           ? "a1cf69ece4dbde9e8bef5b5d80fff4c50c0eaa6e55f43025af4740c7a27c9524" :
    version == "KiLOW_noHe"         ? "5b3de3d3aaad2964f3eeaa8c432bbfca2f2a028c8d62e2f71a64660911a1dc78" :
    version == "KvHIGH_He"          ? "7eb5d3f727c0292ebe91983d6ab8da1f471158cf7e86378f1cde3b487a12f7a8" :
    version == "KvHIGH_KiHIGH_noHe" ? "2037df4277f635d729f1ec681ffcb4b59f818cde0a9ece14ad59c8c3a966988a" :
    version == "KvHIGH_KiLOW_He"    ? "65fc71a09eb019dc3ae54fa069184f3449e5530eb57f68e6e4261786cf74ba5e" :
    version == "KvHIGH_KiLOW_noHe"  ? "bbedae950868639e295056c6b42629593104bf92d55e906b6477d0d77d2b6326" :
    version == "KvHIGH_noHe"        ? "b9a46234e1bcb2da14de6d4fb6e51d7a4ae27b1efbf7df26b47b08718cb30748" :
    OCIM2versionerror(version)
end

# Create registry entry for OCIM in JLD2 format
function register_OCIM2(; version="CTL_He")
    register(
        DataDep(
            "AIBECS-OCIM2_$version",
            """
            References:
            - $(citation())
            """,
            url(version=version),
            sha(version=version),
            fetch_method = fallback_download
        )
    )
    return nothing
end

"""
    load

Returns the grid, the transport matrix, and the He fluxes (in that order).

!!! tip
    To load the default OCIM2 matrix and grid, do
    ```
    julia> grd, T = OCIM2.load()
    ```
    But you can also load the other matrices by specifying which version you want, e.g.,
    ```
    julia> grd, T = OCIM2.load(version="KiHIGH_noHe")
    ```
    See *DeVries and Holzer* (2019) for more details
"""
function load(; version="CTL_He")
    register_OCIM2(version=version)
    jld2_file = @datadep_str string("AIBECS-OCIM2_$version/", "OCIM2_$version.jld2")
    @info """You are about to use the OCIM2_$version model.
          If you use it for research, please cite:

          - $(citation())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "DeVries_Holzer_2019" key.)
          """
    jldopen(jld2_file) do file
        file["grid"], ustrip.(file["T"]), file["He3Flux"], file["He4Flux"]
    end
end

citation() = "DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle‐³He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716"

versions() = ["CTL_He", "CTL_noHe", "KiHIGH_He", "KiHIGH_noHe", "KiLOW_He", "KiLOW_noHe", "KvHIGH_He", "KvHIGH_KiHIGH_noHe", "KvHIGH_KiLOW_He", "KvHIGH_KiLOW_noHe", "KvHIGH_noHe"]

end # end module

export OCIM2

