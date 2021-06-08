module OCIM0
#=
This module serves to load the OCIM0 matrix and grid.
These are loaded from the public, persistant URL in FigShare.
The Julia JLD2 format version of the OCIM0 was created with the code
in the GitHub repository https://github.com/briochemc/OceanCirculations.
=#

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using Downloads
using JLD2                  # For saving circulation as JLD2 format
using CodecZlib             # for JLD2 compression JLD2
using Unitful               # for units
using Reexport
@reexport using OceanGrids            # To store the grid

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Downloads.download(remotepath, localpath)
    return localpath
end

# OCIM0 URLs
function url(; version="")
    url_start = "https://files.figshare.com/"
    return "$(url_start)28336086/OCIM0.1.jld2"
end

# OCIM0 Hashes
sha(; version="") = "4aaaccead0422534ff1adbd208962a8dce7f6cd1af710dcb60faff236863b266"

# Create registry entry for OCIM0 in JLD2 format
function register_OCIM0(; version="")
    register(
        DataDep(
            "AIBECS-OCIM0.1",
            """
            References:
            - $(citation())
            """,
            url(),
            sha(),
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
julia> grd, T = OCIM0.load()
```

See *DeVries and Primeau* (2011) and *Primeau et al.* (2013) for more details.
"""
function load(; version="")
    register_OCIM0()
    jld2_file = @datadep_str string("AIBECS-OCIM0.1/", "OCIM0.1.jld2")
    @info """You are about to use the OCIM0.1 model.
          If you use it for research, please cite:

          $(citation())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "DeVries_Primeau_2011" and "Primeau_etal_2013" keys.)
          """
    jldopen(jld2_file) do file
        file["grid"], ustrip.(file["T"])
    end
end

citation() = """
             - Primeau, F. W., Holzer, M., and DeVries, T. (2013), Southern Ocean nutrient trapping and the efficiency of the biological pump, J. Geophys. Res. Oceans, 118, 2547–2564, doi:10.1002/jgrc.20181.
             - DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, doi:10.1175/JPO-D-10-05011.1
             """

versions() = [""]

end # end module

export OCIM0

