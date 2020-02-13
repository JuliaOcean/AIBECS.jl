module OCIM2
#=
This module serves to load the OCIM2 matrix and grid.
These are loaded from the public, persistant URL in FigShare.
The Julia BSON format version of the OCIM2 was created with the code
in the GitHub repository https://github.com/briochemc/OceanCirculations.
=#

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using BSON                  # For saving circulation as BSON format
using Unitful, UnitfulAstro # for units
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
function register_OCIM2()
    register(
        DataDep(
            "AIBECS_OCIM2",
            """
            References:
            - DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle‐³He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716 
            """,
            "https://files.figshare.com/21689934/OCIM2.bson",
            # use `sha2_256` (no quotes) to figure out the hash on first try
            "77b018c75006409ecd0e44a7371033e9b554ba0cd50506014bc8635f12be5042",
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
    print("Loading OCIM2")
    register_OCIM2()
    bson_file = @datadep_str string("AIBECS_OCIM2/", "OCIM2.bson")
    BSON.@load bson_file T grid
    println(" ✔")
    @info """You are about to use the OCIM2 (CTL with He) model.
          If you use it for research, please cite:

          - DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle‐³He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716 

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "DeVries_Holzer_2019" key.)
          """
    return grid, ustrip.(T)
end

end # end module

export OCIM2

