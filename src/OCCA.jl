"""
This `OCCA` module is used to load the OCCA matrix and grid for use in AIBECS.

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

# OCCA URLs
function url(; version="OCCA")
    url_start = "https://files.figshare.com/"
    version == "OCCA"             ? "$(url_start)22823813/OCCA.bson" :
    OCCAversionerror(version)
end

OCCAversionerror(version) = error("""`$version` is not a valid OCCA circulation name.

                                   Valid versions are `OCCA`.

                                   See *Forget* (2010) for more details on OCCA.""")

# OCCA Hashes
function sha(; version="OCCA")
    version == "OCCA" ? sha2_256 :
    OCCAversionerror(version)
end

# Create registry entry for OCCA in BSON format
function register_OCCA(; version="OCCA")
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

Returns the grid and the transport matrix.

!!! tip
    To load the OCCA matrix and grid, do
    ```
    julia> grd, T = OCCA.load()
    ```
    See *Forget* (2010) for more details
"""
function load(; version="OCCA")
    register_OCCA(version=version)
    bson_file = @datadep_str string("AIBECS-$version/", "$version.bson")
    BSON.@load bson_file grid T
    @info """You are about to use the $version model.
          If you use it for research, please cite:

          - $(citation())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "Forget_2010" key.)
          """
    return grid, ustrip.(T)
end

citation() = "Forget, G., 2010: Mapping Ocean Observations in a Dynamical Framework: A 2004–06 Ocean Atlas. J. Phys. Oceanogr., 40, 1201–1221, https://doi.org/10.1175/2009JPO4043.1"

end # end module

export OCCA

