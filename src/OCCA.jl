"""
The `OCCA` module is used to load the OCCA matrix and grid for use in AIBECS.

!!! tip
    To load the OCCA matrix and grid, do
    ```
    julia> grd, T = OCCA.load()
    ```
    See *Forget* (2010) for more details

!!! note
    The files, that are downloaded from a public and persistent URL on Zenodo,
    were created with the code available at https://github.com/briochemc/OceanCirculations.

!!! warning "Mass conservation"
    The raw OCCA matrix is not mass-conserving, which can cause issues.
    To avoid these issues, AIBECS's `load()` call will enforce mass conservation
    by default. This is done by subtracting `Diagonal((Tᵀ v) ./ v)` from
    `T`. Pass `conservemass = false` to recover the raw matrix; in that case
    tracer budgets will drift unless the geological restoring is modeled
    explicitly.
"""
module OCCA

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using Downloads
using Unitful               # for units
using Reexport
using MD5                   # for hash checking (MD5 is what Zenodo exposes)
@reexport using OceanGrids  # To store the grid

# OCCA URLs
const OCCA_URL = "https://zenodo.org/records/20079672/files/OCCA.jld2"

# OCCA Hashes
const OCCA_MD5 = "bde4109fb20ee876cfdfc68d2a10ba76"

const CITATION = "Forget, G., 2010: Mapping Ocean Observations in a Dynamical Framework: A 2004–06 Ocean Atlas. J. Phys. Oceanogr., 40, 1201–1221, https://doi.org/10.1175/2009JPO4043.1"

# Create registry entry for OCCA in JLD2 format
function register_OCCA()
    register(
        DataDep(
            "AIBECS-OCCA",
            """
            References:
            - $CITATION
            """,
            OCCA_URL,
            (md5, OCCA_MD5)
        )
    )
    return nothing
end

invalidate_stale_cache() =
    parentmodule(@__MODULE__)._invalidate_stale_cache("AIBECS-OCCA", "OCCA.jld2", OCCA_MD5)

"""
    grd, T = load()

Returns the grid and the transport matrix.

Requires `using JLD2` so that the `AIBECSJLD2Ext` extension is activated.

See *Forget* (2010) for more details.
"""
function load(args...; kwargs...)
    error("AIBECS.OCCA.load requires `using JLD2`. Add it to your environment, then retry.")
end

end # end module
