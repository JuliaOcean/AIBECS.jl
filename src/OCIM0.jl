"""
The `OCIM0` module is used to load the OCIM v0 matrix and grid for use in AIBECS.

!!! tip
    To load the OCIM v0 matrix and grid, do
    ```
    julia> grd, T = OCIM0.load()
    ```
    See *DeVries and Primeau* (2011) and *Primeau et al.* (2013) for more details.

!!! note
    The files, that are downloaded from a public and persistent URL on Zenodo,
    were created with the code available at https://github.com/briochemc/OceanCirculations.
"""
module OCIM0

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using Downloads
using Unitful               # for units
using Reexport
using MD5                   # for hash checking (MD5 is what Zenodo exposes)
@reexport using OceanGrids  # To store the grid

# OCIM0 URLs
const OCIM0_URL = "https://zenodo.org/records/19944573/files/OCIM0.1.jld2"

# OCIM0 Hashes
const OCIM0_MD5 = "458a463a9f1fdc223a6cfb025daa2d47"

const CITATION = """
- Primeau, F. W., Holzer, M., and DeVries, T. (2013), Southern Ocean nutrient trapping and the efficiency of the biological pump, J. Geophys. Res. Oceans, 118, 2547–2564, doi:10.1002/jgrc.20181.
- DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, doi:10.1175/JPO-D-10-05011.1
"""

# Create registry entry for OCIM0 in JLD2 format
function register_OCIM0()
    register(
        DataDep(
            "AIBECS-OCIM0.1",
            """
            References:
            - $CITATION
            """,
            OCIM0_URL,
            (md5, OCIM0_MD5),
        )
    )
    return nothing
end

"""
    grd, T = load()

Returns the grid and the transport matrix.

Requires `using JLD2` so that the `AIBECSJLD2Ext` extension is activated.

See *DeVries and Primeau* (2011) and *Primeau et al.* (2013) for more details.
"""
function load(args...; kwargs...)
    error("AIBECS.OCIM0.load requires `using JLD2`. Add it to your environment, then retry.")
end


end # end module

export OCIM0

