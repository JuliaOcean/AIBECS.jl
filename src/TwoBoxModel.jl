module TwoBoxModel
#=
This module serves to load the 2-box model matrix and grid
from the Sarmiento and Gruber book (2006) book.

┌────────────────────────────────┐
│                        Surface │
├────────────────────────────────┤
│                                │
│                                │
│                                │
│                                │
│                           Deep │
└────────────────────────────────┘

=#

using LinearAlgebra, SparseArrays
using Unitful, UnitfulAstro # for units
using Reexport
@reexport using OceanGrids            # To store the grid
using ..CirculationGeneration
CG = CirculationGeneration

build_wet3D() = trues(1, 1, 2)

function build_grid()
    wet3D = build_wet3D()
    vdeep = 1.26e18u"m^3"     # Sarmiento and Gruber, 2006
    vtot = 1.292e18u"m^3"     # Archer et al., 2000
    area_ocean = 349e12u"m^2" # Archer et al., 2000
    # Infer other values
    depth_surface = (vtot-vdeep) / area_ocean
    depth_ocean = vtot / area_ocean
    # Build DY3d
    elon = [0,360] * u"°"
    elat = [-90, 90] * u"°"
    edepth = [0u"m", depth_surface, depth_ocean]
    return OceanGrid(elon, elat, edepth, wet3D)
end

function build_T(grid)
    # From Archer et al. [2000]
    v3D = array_of_volumes(grid)
    nb = length(v3D)
    # ν = 38 Sv (Sarmiento and Gruber, 2008) (1 Sv = 10⁶ m³/s)
    T  = CG.flux_divergence_operator_from_advection(38e6u"m^3/s", [1, 2], v3D, nb)
    return T
end



"""
    load

Returns `grd` and `T` (in that order).
"""
function load()
    print("Creating 2-box model")
    grid = build_grid()
    T = build_T(grid)
    println(" ✔")
    @info """
            You are about to use the 2-box model of Sarmiento and Gruber [2006].
            If you use it for research, please cite:

            - Sarmiento, Jorge L., and Nicolas Gruber. Ocean biogeochemical dynamics. Princeton University Press, 2006.

            You can find the corresponding BibTeX entry in the CITATION.bib file
            at the root of the AIBECS.jl package repository.
            (Look for the "Sarmiento_Gruber_2006" key.)
            """
    return grid, T
end

end

export TwoBoxModel

