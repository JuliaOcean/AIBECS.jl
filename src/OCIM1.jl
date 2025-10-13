"""
The `OCIM1` module is used to load the OCIM v1 matrix and grid for use in AIBECS.

!!! tip
    To load the OCIM v1 matrix and grid, do
    ```
    julia> grd, T = OCIM1.load()
    ```
    See *DeVries and Primeau* (2011) and *Primeau et al.* (2013) for more details.

!!! note
    The files, that are downloaded from a public and persistant URL in FigShare,
    were created with the code available at https://github.com/briochemc/OceanCirculations.
"""
module OCIM1

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using Downloads
using JLD2                  # For saving circulation as JLD2 format
using CodecZlib             # for JLD2 compression JLD2
using Unitful               # for units
using Reexport
using MD5                   # for hash checking (MD5 is what is used in FigShare)
@reexport using OceanGrids  # To store the grid

const VERSIONS = ["CTL"] # TODO add different OCIM1 versions

# URLs
const OCIM1_URL = "https://ndownloader.figshare.com/files/28335882"

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
    load

Returns the grid and the transport matrix.

!!! tip
    To load the OCIM1 matrix and grid, do
    ```
    julia> grd, T = OCIM1.load()
    ```
    See *DeVries and Primeau* (2011) and *DeVries* (2014) for more details.
"""
function load(; version=VERSIONS[1])
    register_OCIM1(; version)
    jld2_file = @datadep_str string("AIBECS-OCIM1_$version/", "OCIM1_$version.jld2")
    @info """You are about to use the OCIM1_$version model.
          If you use it for research, please cite:

          $CITATION

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "DeVries_Primeau_2011" and "DeVries_2014" keys.)
          """
    jldopen(jld2_file) do file
        file["grid"], ustrip.(file["T"])
    end
end

end # end module

export OCIM1

