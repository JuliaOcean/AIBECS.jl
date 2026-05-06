"""
The `OCIM2_48L` module is used to load the OCIM2-48L matrices and grid for use in AIBECS.

!!! tip
    To load the default OCIM2_48L matrix and grid, do
    ```
    julia> grd, T = OCIM2_48L.load()
    ```
    But you can also load the other matrices by specifying which version you want, e.g.,
    ```
    julia> grd, T = OCIM2_48L.load(version="OCIM2_48L_KiHIGH_noHe")
    ```
    See *DeVries and Holzer* (2019) for more details


!!! note
    The files, that are downloaded from a public and persistent URL on Zenodo,
    were created with the code available at https://github.com/briochemc/OceanCirculations.
"""
module OCIM2_48L

using SparseArrays          # for sparse matrix
using DataDeps              # for automated download/storage of data
using Downloads
using Unitful               # for units
using Unitful: s, yr, m, km, °
using MD5
using Reexport
@reexport using OceanGrids            # To store the grid

# OCIM2_48L URL
const URL = "https://zenodo.org/records/19944665/files/OCIM2_48L_base.tar.gz"

# OCIM2_48L Hash
const MD5 = "7938f20f06eff072b2791b571d2bb9d7"

OCIM2versionerror(version) = error(
    """`$version` is not a valid OCIM2 version name.

    Valid versions are `CTL_He`, `CTL_noHe`, `KiHIGH_He`, `KiHIGH_noHe`, `KiLOW_He`, `KiLOW_noHe`, `KvHIGH_He`, `KvHIGH_KiHIGH_noHe`, `KvHIGH_KiLOW_He`, `KvHIGH_KiLOW_noHe`, and `KvHIGH_noHe`.

    See *DeVries and Holzer* (2019) for more details on OCIM2 configurations."""
)


# Create registry entry for OCIM in JLD2 format
function register_OCIM2_48L()
    register(
        DataDep(
            "AIBECS-OCIM2_48L",
            """
            References:
            - $CITATION
            """,
            URL,
            (md5, MD5),
            post_fetch_method = unpack
        )
    )
    return nothing
end

"""
    grd, T = load()

Returns the grid and the transport matrix of the OCIM2_48L model.

See *DeVries and Holzer* (2019) and *Holzer et al.* (2021) for more details.

Requires `using MAT, NCDatasets` so that the `AIBECSOCIM2_48LExt` extension is activated.
"""
function load(args...; kwargs...)
    error(
        "AIBECS.OCIM2_48L.load requires `using MAT, NCDatasets`. " *
            "Add them to your environment, then retry."
    )
end

const CITATION = """
- Holzer, M., DeVries, T. & de Lavergne, C. Diffusion controls the ventilation of a Pacific Shadow Zone above abyssal overturning. Nat Commun 12, 4348 (2021). https://doi.org/10.1038/s41467-021-24648-x
- DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle-3He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716
"""

end # end module
