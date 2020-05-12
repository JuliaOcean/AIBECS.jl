module OCIM0
#=
This module serves to load the OCIM0 matrix and grid.
These are loaded from the public, persistant URL in FigShare.
The Julia BSON format version of the OCIM0 was created with the code
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

# Create registry entry for OCIM in JLD2 format
function register_OCIM0()
    register(
        DataDep(
            "AIBECS_OCIM0.1",
            """
            References:
            - Primeau, F. W., Holzer, M., and DeVries, T. (2013), Southern Ocean nutrient trapping and the efficiency of the biological pump, J. Geophys. Res. Oceans, 118, 2547–2564, doi:10.1002/jgrc.20181.
            - DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, https://doi.org/10.1175/JPO-D-10-05011.1
            """,
            "https://files.figshare.com/18789281/OCIM0.1.bson",
            # use `sha2_256` (no quotes) to figure out the hash on first try
            "527f02545ecafb59f78eeaa616c012274608971fba887eea1384aa1d9b0d40a9",
            fetch_method = fallback_download
        )
    )
    return nothing
end

"""
    load

Returns `grd` and `T` (in that order) from FigShare repository.
"""
function load()
    print("Loading OCIM0.1")
    register_OCIM0()
    bson_file = @datadep_str string("AIBECS_OCIM0.1/", "OCIM0.1.bson")
    BSON.@load bson_file T grid
    println(" ✔")
    @info """You are about to use OCIM0.1 model.
          If you use it for research, please cite:

          - Primeau, F. W., Holzer, M., and DeVries, T. (2013), Southern Ocean nutrient trapping and the efficiency of the biological pump, J. Geophys. Res. Oceans, 118, 2547–2564, doi:10.1002/jgrc.20181.
          - DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, doi:10.1175/JPO-D-10-05011.1

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "DeVries_Primeau_2011" and "Primeau_etal_2013" keys.)
          """
    return grid, ustrip.(T)
end

end # end module

export OCIM0

