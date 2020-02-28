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
using BSON                  # For saving circulation as BSON format
using Unitful, UnitfulAstro # for units
using Reexport
@reexport using OceanGrids            # To store the grid

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Base.download(remotepath, localpath)
    return localpath
end

# OCIM2 URLs
function url(; version="OCIM2_CTL_He")
    url_start = "https://files.figshare.com/"
    version == "OCIM2_CTL_He"             ? "$(url_start)21855225/OCIM2_CTL_He.bson" :
    version == "OCIM2_CTL_noHe"           ? "$(url_start)21843960/OCIM2_CTL_noHe.bson" :
    version == "OCIM2_KiHIGH_He"          ? "$(url_start)21843963/OCIM2_KiHIGH_He.bson" :
    version == "OCIM2_KiHIGH_noHe"        ? "$(url_start)21843966/OCIM2_KiHIGH_noHe.bson" :
    version == "OCIM2_KiLOW_He"           ? "$(url_start)21843969/OCIM2_KiLOW_He.bson" :
    version == "OCIM2_KiLOW_noHe"         ? "$(url_start)21843972/OCIM2_KiLOW_noHe.bson" :
    version == "OCIM2_KvHIGH_He"          ? "$(url_start)21843984/OCIM2_KvHIGH_He.bson" :
    version == "OCIM2_KvHIGH_KiHIGH_noHe" ? "$(url_start)21843987/OCIM2_KvHIGH_KiHIGH_noHe.bson" :
    version == "OCIM2_KvHIGH_KiLOW_He"    ? "$(url_start)21843990/OCIM2_KvHIGH_KiLOW_He.bson" :
    version == "OCIM2_KvHIGH_KiLOW_noHe"  ? "$(url_start)21843993/OCIM2_KvHIGH_KiLOW_noHe.bson" :
    version == "OCIM2_KvHIGH_noHe"        ? "$(url_start)21843996/OCIM2_KvHIGH_noHe.bson" :
    OCIM2versionerror(version)
end

OCIM2versionerror(version) = error("""`$version` is not a valid OCIM2 circulation name.

                                   Valid versions are `OCIM2_CTL_He`, `OCIM2_CTL_noHe`, `OCIM2_KiHIGH_He`, `OCIM2_KiHIGH_noHe`, `OCIM2_KiLOW_He`, `OCIM2_KiLOW_noHe`, `OCIM2_KvHIGH_He`, `OCIM2_KvHIGH_KiHIGH_noHe`, `OCIM2_KvHIGH_KiLOW_He`, `OCIM2_KvHIGH_KiLOW_noHe`, and `OCIM2_KvHIGH_noHe`.

                                   See *DeVries and Holzer* (2019) for more details on OCIM2 configurations.""")

# OCIM2 Hashes
function sha(; version="OCIM2_CTL_He")
    version == "OCIM2_CTL_He"             ? "bd3bb098543740872025325b6b16a78bd88091be2778a3d394eefc69d1d70889" :
    version == "OCIM2_CTL_noHe"           ? "bdfbdc67b224ea8a588450454feb1f440cfb3745c191ff3929ac66897a036e02" :
    version == "OCIM2_KiHIGH_He"          ? "5ad8f8d2b1fd3c65bc89999e7a0bb09e178939d91983d8611bdf9d26a7497ef7" :
    version == "OCIM2_KiHIGH_noHe"        ? "9ec3f7237e541185bbf577bc7bee8c6308f2e4a2e0457a7460e73a962902cd8c" :
    version == "OCIM2_KiLOW_He"           ? "b26e35eca25b48a211d94b9a4cf0f260b9ae5a6a0acad673b1060ff6e7522e69" :
    version == "OCIM2_KiLOW_noHe"         ? "a732a382365235446091f77c2cf2961af641b10662a653edec05295349151e8f" :
    version == "OCIM2_KvHIGH_He"          ? "92ae46c53ae632ada056c6333fba885a758d72be40eefedf87c08f60d18a8c60" :
    version == "OCIM2_KvHIGH_KiHIGH_noHe" ? "8ee547b2fd4ddc0b61796f0c4bc0e57622615c1e492939377b44a899f86a3122" :
    version == "OCIM2_KvHIGH_KiLOW_He"    ? "ce362eef6e5265d411d79fdfe3abba40498f399d6f490ff3dd51fb2b163595bf" :
    version == "OCIM2_KvHIGH_KiLOW_noHe"  ? "22c3de72f6d630901d929c769d1c88c09ba15264f489c34c1b0d1069fc2538c5" :
    version == "OCIM2_KvHIGH_noHe"        ? "419cc207a0e104427677b816c645e8fdbca83d1201204dd4dd10a95200fbc070" :
    OCIM2versionerror(version)
end

# Create registry entry for OCIM in JLD2 format
function register_OCIM2(; version="OCIM2_CTL_He")
    register(
        DataDep(
            "AIBECS-$version",
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
    julia> grd, T = OCIM2.load(version="OCIM2_KiHIGH_noHe")
    ```
    See *DeVries and Holzer* (2019) for more details
"""
function load(; version="OCIM2_CTL_He")
    register_OCIM2(version=version)
    bson_file = @datadep_str string("AIBECS-$version/", "$version.bson")
    BSON.@load bson_file grid T He3Flux He4Flux
    @info """You are about to use the $version model.
          If you use it for research, please cite:

          - $(citation())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "DeVries_Holzer_2019" key.)
          """
    return grid, ustrip.(T), He3Flux, He4Flux
end

citation() = "DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle‐³He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716"

end # end module

export OCIM2

