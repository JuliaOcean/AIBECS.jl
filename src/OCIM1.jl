"""
The `OCIM1` module is used to load the OCIM v1 matrix and grid for use in AIBECS.

!!! tip
    To load the OCIM v1 matrix and grid, do
    ```
    julia> grd, T = OCIM1.load()
    ```
    See *DeVries and Primeau* (2011) and *Primeau et al.* (2013) for more details.

!!! note
    The files, that are downloaded from a public and persistent URL on Zenodo,
    were created with the code available at https://github.com/briochemc/OceanCirculations.
"""
module OCIM1

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using Downloads
using Unitful               # for units
using Reexport
using MD5                   # for hash checking (MD5 is what Zenodo exposes)
@reexport using OceanGrids  # To store the grid

const VERSIONS = ["CTL"] # TODO add different OCIM1 versions

# URLs
const OCIM1_URL = "https://zenodo.org/records/19944597/files/OCIM1_CTL.jld2"

# OCIM1 Hashes
const OCIM1_MD5 = "eaa57b42e7edec0fe965575e9938c66d"

const CITATION = """
- DeVries, T., 2014: The oceanic anthropogenic CO2 sink: Storage, air‐sea fluxes, and transports over the industrial era, Global Biogeochem. Cycles, 28, 631–647, doi:10.1002/2013GB004739.
- DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, doi:10.1175/JPO-D-10-05011.1
"""

# Create registry entry for OCIM1 in JLD2 format
function register_OCIM1(; version=VERSIONS[1])
    register(
        DataDep(
            "AIBECS-OCIM1_$version",
            """
            References:
            $CITATION
            """,
            OCIM1_URL,
            (md5, OCIM1_MD5),
        )
    )
    return nothing
end

"""
    grd, T = load(; version=VERSIONS[1])

Returns the grid and the transport matrix.

Requires `using JLD2` so that the `AIBECSJLD2Ext` extension is activated.

See *DeVries and Primeau* (2011) and *DeVries* (2014) for more details.
"""
function load(args...; kwargs...)
    error("AIBECS.OCIM1.load requires `using JLD2`. Add it to your environment, then retry.")
end

end # end module

export OCIM1

