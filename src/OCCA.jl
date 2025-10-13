"""
The `OCCA` module is used to load the OCCA matrix and grid for use in AIBECS.

!!! tip
    To load the OCCA matrix and grid, do
    ```
    julia> grd, T = OCCA.load()
    ```
    See *Forget* (2010) for more details

!!! note
    The files, that are downloaded from a public and persistant URL in FigShare,
    were created with the code available at https://github.com/briochemc/OceanCirculations.
"""
module OCCA

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using Downloads
using JLD2                  # For saving circulation as JLD2 format
using CodecZlib             # for JLD2 compression JLD2
using Unitful               # for units
using Reexport
using MD5                   # for hash checking (MD5 is what is used in FigShare)
@reexport using OceanGrids  # To store the grid

# OCCA URLs
const OCCA_URL = "https://ndownloader.figshare.com/files/28336173"

# OCCA Hashes
const OCCA_MD5 = "f2a2a2d295f85771c2302e6c0eb35f4c"

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
            (md5, OCCA_MD5);
        )
    )
    return nothing
end

"""
    load

Returns the grid and the transport matrix.

!!! tip
    To load the OCCA matrix and grid, do
    ```
    julia> grd, T = OCCA.load()
    ```
    See *Forget* (2010) for more details
"""
function load()
    register_OCCA()
    jld2_file = @datadep_str string("AIBECS-OCCA/", "OCCA.jld2")
    @info """You are about to use the OCCA model.
          If you use it for research, please cite:

          - $CITATION

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "Forget_2010" key.)
          """
    jldopen(jld2_file) do file
        file["grid"], ustrip.(file["T"])
    end
end

end # end module

export OCCA
