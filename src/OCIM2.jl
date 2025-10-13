"""
The `OCIM2` module is used to load the OCIM2 matrices and grid for use in AIBECS.

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
using MD5                   # for hash checking (MD5 is what is used in FigShare)
@reexport using OceanGrids  # To store the grid

const VERSIONS =[
    "CTL_He",
    "CTL_noHe",
    "KiHIGH_He",
    "KiHIGH_noHe",
    "KiLOW_He",
    "KiLOW_noHe",
    "KvHIGH_He",
    "KvHIGH_KiHIGH_noHe",
    "KvHIGH_KiLOW_He",
    "KvHIGH_KiLOW_noHe",
    "KvHIGH_noHe"
]

# URLs
const OCIM1_URLs = Dict(
    "CTL_He" => "https://ndownloader.figshare.com/files/28336284",
    "CTL_noHe" => "https://ndownloader.figshare.com/files/28336299",
    "KiHIGH_He" => "https://ndownloader.figshare.com/files/28336302",
    "KiHIGH_noHe" => "https://ndownloader.figshare.com/files/28336311",
    "KiLOW_He" => "https://ndownloader.figshare.com/files/28336317",
    "KiLOW_noHe" => "https://ndownloader.figshare.com/files/28336323",
    "KvHIGH_He" => "https://ndownloader.figshare.com/files/28336326",
    "KvHIGH_KiHIGH_noHe" => "https://ndownloader.figshare.com/files/28336329",
    "KvHIGH_KiLOW_He" => "https://ndownloader.figshare.com/files/28336332",
    "KvHIGH_KiLOW_noHe" => "https://ndownloader.figshare.com/files/28336341",
    "KvHIGH_noHe" => "https://ndownloader.figshare.com/files/28336353"
)

# OCIM1 Hashes
const OCIM2_MD5 = Dict(
    "CTL_He" => "da15192381ef04e9f7cce7c886eb4833",
    "CTL_noHe" => "abede7e16b75cb3f9b7d25367772455e",
    "KiHIGH_He" => "3c8312f9aeb9bb6cb6ac32a8d84fb901",
    "KiHIGH_noHe" => "781b8c058eb74a81d10befa83faf0d08",
    "KiLOW_He" => "ee6c6f14c6bd90c1980daf8a220f0551",
    "KiLOW_noHe" => "47eabd73588739344c81f72d39367f27",
    "KvHIGH_He" => "ecd18348d062fb1c1e1518244b2b8275",
    "KvHIGH_KiHIGH_noHe" => "c71514b13802c0e497f033612f586f5a",
    "KvHIGH_KiLOW_He" => "99edb2894a2ee274e68f3309da1ed200",
    "KvHIGH_KiLOW_noHe" => "ff8fb1a0add22bde6f3fed393fa7336f",
    "KvHIGH_noHe" => "0c0402b3796c3a7d14c0bdc2db7f0aa7"
)

const CITATION = """
DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle‐³He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716"
"""


# Create registry entry for OCIM2 in JLD2 format
function register_OCIM2(; version=VERSIONS[1])
    register(
        DataDep(
            "AIBECS-OCIM2_$version",
            """
            References:
            $CITATION
            """,
            OCIM2_URLs[version],
            (md5, OCIM1_MD5[version]),
        )
    )
    return nothing
end

"""
    load

Returns the grid, the transport matrix, and the He fluxes.

!!! tip
    To load the default OCIM2 matrix and grid, do
    ```
    julia> grd, T = OCIM2.load()
    ```
    But you can also load the other matrices by specifying which version you want, e.g.,
    ```
    julia> grd, T = OCIM2.load(version="KiHIGH_noHe")
    ```

    Add the `HeFluxes=true` keyword argument if you want the OCIM-produced He fields
    with ³He and ⁴He as 3rd and 4th arguments.
    (3rd and 4th output are returned as `nothing` for "noHe" versions).

    See *DeVries and Holzer* (2019) for more details.

    To see all available versions, do
    ```
    julia> OCIM2.VERSIONS
    ```
"""
function load(; version=VERSIONS[1], HeFluxes=false)
    register_OCIM2(; version)
    jld2_file = @datadep_str string("AIBECS-OCIM2_$version/", "OCIM2_$version.jld2")
    @info """You are about to use the OCIM2_$version model.
          If you use it for research, please cite:

          $CITATION

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "DeVries_Holzer_2019" key.)
          """
    if HeFluxes
        jldopen(jld2_file) do file
            file["grid"], ustrip.(file["T"]), file["He3Flux"], file["He4Flux"]
        end
    else
        jldopen(jld2_file) do file
            file["grid"], ustrip.(file["T"])
        end
    end
end

end # end module

export OCIM2

