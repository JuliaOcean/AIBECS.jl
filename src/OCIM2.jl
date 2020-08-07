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
using Unitful               # for units
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
function url(; version="CTL_He")
    url_start = "https://files.figshare.com/"
    version == "CTL_He"             ? "$(url_start)22475534/OCIM2_CTL_He.bson" :
    version == "CTL_noHe"           ? "$(url_start)23442083/OCIM2_CTL_noHe.bson" :
    version == "KiHIGH_He"          ? "$(url_start)23419721/OCIM2_KiHIGH_He.bson" :
    version == "KiHIGH_noHe"        ? "$(url_start)23418215/OCIM2_KiHIGH_noHe.bson" :
    version == "KiLOW_He"           ? "$(url_start)23416925/OCIM2_KiLOW_He.bson" :
    version == "KiLOW_noHe"         ? "$(url_start)23417861/OCIM2_KiLOW_noHe.bson" :
    version == "KvHIGH_He"          ? "$(url_start)23416583/OCIM2_KvHIGH_He.bson" :
    version == "KvHIGH_KiHIGH_noHe" ? "$(url_start)23417639/OCIM2_KvHIGH_KiHIGH_noHe.bson" :
    version == "KvHIGH_KiLOW_He"    ? "$(url_start)23416181/OCIM2_KvHIGH_KiLOW_He.bson" :
    version == "KvHIGH_KiLOW_noHe"  ? "$(url_start)23416055/OCIM2_KvHIGH_KiLOW_noHe.bson" :
    version == "KvHIGH_noHe"        ? "$(url_start)23417396/OCIM2_KvHIGH_noHe.bson" :
    OCIM2versionerror(version)
end

OCIM2versionerror(version) = error("""`$version` is not a valid OCIM2 version name.

                                   Valid versions are `CTL_He`, `CTL_noHe`, `KiHIGH_He`, `KiHIGH_noHe`, `KiLOW_He`, `KiLOW_noHe`, `KvHIGH_He`, `KvHIGH_KiHIGH_noHe`, `KvHIGH_KiLOW_He`, `KvHIGH_KiLOW_noHe`, and `KvHIGH_noHe`.

                                   See *DeVries and Holzer* (2019) for more details on OCIM2 configurations.""")

# OCIM2 Hashes
function sha(; version="CTL_He")
    version == "CTL_He"             ? "2786a56e89315feb44102629470ec565bd67ee765cb0a49364a51e6206a15490" :
    version == "CTL_noHe"           ? "bdfbdc67b224ea8a588450454feb1f440cfb3745c191ff3929ac66897a036e02" :
    version == "KiHIGH_He"          ? "bf89d657b690bbfc8055c2182c04d50f5a408b7c610b4dc8daad52fba2482bf4" :
    version == "KiHIGH_noHe"        ? "9ec3f7237e541185bbf577bc7bee8c6308f2e4a2e0457a7460e73a962902cd8c" :
    version == "KiLOW_He"           ? "2a7094a5fc0bf5aa6d5d4e0f7ad4ed89901a7431a8d039df74835d84d4259a91" :
    version == "KiLOW_noHe"         ? "a732a382365235446091f77c2cf2961af641b10662a653edec05295349151e8f" :
    version == "KvHIGH_He"          ? "bf89d657b690bbfc8055c2182c04d50f5a408b7c610b4dc8daad52fba2482bf4" :
    version == "KvHIGH_KiHIGH_noHe" ? "8ee547b2fd4ddc0b61796f0c4bc0e57622615c1e492939377b44a899f86a3122" :
    version == "KvHIGH_KiLOW_He"    ? "c207e88e88a155972f51a82e1ec216efb9bd28073464074a7b023460ff9805b6" :
    version == "KvHIGH_KiLOW_noHe"  ? "22c3de72f6d630901d929c769d1c88c09ba15264f489c34c1b0d1069fc2538c5" :
    version == "KvHIGH_noHe"        ? "419cc207a0e104427677b816c645e8fdbca83d1201204dd4dd10a95200fbc070" :
    OCIM2versionerror(version)
end

# Create registry entry for OCIM in BSON format
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
    bson_file = @datadep_str string("AIBECS-OCIM2_$version/", "OCIM2_$version.bson")
    BSON.@load bson_file grid T He3Flux He4Flux
    @info """You are about to use the OCIM2_$version model.
          If you use it for research, please cite:

          - $(citation())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "DeVries_Holzer_2019" key.)
          """
    return grid, ustrip.(T), He3Flux, He4Flux
end

citation() = "DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle‐³He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716"

versions() = ["CTL_He", "CTL_noHe", "KiHIGH_He", "KiHIGH_noHe", "KiLOW_He", "KiLOW_noHe", "KvHIGH_He", "KvHIGH_KiHIGH_noHe", "KvHIGH_KiLOW_He", "KvHIGH_KiLOW_noHe", "KvHIGH_noHe"]

end # end module

export OCIM2

