module OCIM1
#=
This module serves to load the OCIM1 matrix and grid.
These are loaded from the public, persistant URL in FigShare.
The Julia BSON format version of the OCIM1 was created with the code
in the GitHub repository https://github.com/briochemc/OceanCirculations.
=#

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

# OCIM1 URLs
function url(; version="CTL") # TODO add other OCIM1 circulations
    url_start = "https://files.figshare.com/"
    version == "CTL" ? "$(url_start)23420027/OCIM1_CTL.bson" :
    OCIM1versionerror(version)
end

OCIM1versionerror(version) = error("""`$version` is not a valid OCIM1 version name (from within AIBECS).

                                   Only valid version is `CTL`.

                                   Start an issue on GitHub if you want me to add other OCIM1 circulations.""")

# OCIM1 Hashes
function sha(; version="CTL")
    version == "CTL" ? "28d4290c04844180191e71d5e1bc832ff43c1127dd621b719c6622a2288d5e24" :
    OCIM1versionerror(version)
end

# Create registry entry for OCIM1 in BSON format
function register_OCIM1(; version="CTL")
    register(
        DataDep(
            "AIBECS-OCIM1_$version",
            """
            References:
            $(citation())
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

Returns the grid and the transport matrix (in that order).

### Usage

```
julia> grd, T = OCIM1.load()
```

See *DeVries and Primeau* (2011) and *DeVries* (2014) for more details.
"""
function load(; version="CTL")
    register_OCIM1(version=version)
    bson_file = @datadep_str string("AIBECS-OCIM1_$version/", "OCIM1_$version.bson")
    BSON.@load bson_file grid T
    @info """You are about to use the OCIM1_$version model.
          If you use it for research, please cite:

          $(citation())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "DeVries_Primeau_2011" and "DeVries_2014" keys.)
          """
    return grid, ustrip.(T)
end

citation() = """
             - DeVries, T., 2014: The oceanic anthropogenic CO2 sink: Storage, air‐sea fluxes, and transports over the industrial era, Global Biogeochem. Cycles, 28, 631–647, doi:10.1002/2013GB004739.
             - DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, doi:10.1175/JPO-D-10-05011.1
             """

versions() = ["CTL"]

end # end module

export OCIM1

