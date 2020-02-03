module Archer_etal_2000
#=
This module serves to load the 6-box model matrix and grid.
The 6-box model is used because it makes it easier to build
operators for sinking particles in a regular grid,
and because the code to bin WOA data requires a regular grid.
These are the indices of the 6-box model:

┌───────┬────────────────────────┐
│ 1(H)  │         2(L)           │
├───────┼────────────────────────┤
│ 3(H)  │         4(D)           │
├───────┼────────────────────────┤
│       │                        │
│ 5(D)  │         6(D)           │
│       │                        │
└───────┴────────────────────────┘

The boxes are grouped as indicated:
- 1 + 3 represents the high-latitudes (surface) box
- 2 represents the low-latitudes (surface) box
- 4 + 5 + 6 represents the deep box
=#

using LinearAlgebra, SparseArrays
using Unitful, UnitfulAstro # for units
using Reexport
@reexport using OceanGrids            # To store the grid
using ..CirculationGeneration
CG = CirculationGeneration

build_wet3D() = trues(2, 1, 3)

function build_grid()
    wet3D = build_wet3D()
    # From Archer et al. [2000]
    vtot = 1.292e18u"m^3"
    area_ocean = 349e12u"m^2"
    depth_low = 100.0u"m"
    depth_high = 250.0u"m"
    fraction_area_high = 0.15 #
    # Infer other values
    depth_ocean = vtot / area_ocean
    # Build DY3d
    elon = [0,360] * u"°"
    elat = [-90, -90+fraction_area_high*180, 90] * u"°"
    edepth = [0u"m", depth_low, depth_high, depth_ocean]
    return OceanGrid(elon, elat, edepth, wet3D)
end

function build_T(grid)
    # From Archer et al. [2000]
    v3D = array_of_volumes(grid)
    nb = length(v3D)

    # Mixing terms
    # between top boxes (1 and 2) (1 Sv = 10⁶ m³/s)
    T  = CG.flux_divergence_operator_from_advection(10e6u"m^3/s", [1, 2], v3D, nb)
    # between high-lat and deep (1 and 3)
    T += CG.flux_divergence_operator_from_advection(53e6u"m^3/s", [2, 5], v3D, nb)
    # between low-lat and deep (2 and 4)
    T += CG.flux_divergence_operator_from_advection( 1e6u"m^3/s", [2, 4], v3D, nb)
    # (trick) fast mixing between boxes of the same 3-box-model box
    T += CG.flux_divergence_operator_from_advection(1e10u"m^3/s", [1, 3], v3D, nb)
    T += CG.flux_divergence_operator_from_advection(1e10u"m^3/s", [4, 5], v3D, nb)
    T += CG.flux_divergence_operator_from_advection(1e10u"m^3/s", [5, 6], v3D, nb)
    T += CG.flux_divergence_operator_from_advection(1e10u"m^3/s", [4, 6], v3D, nb)

    # Overturning circulation
    T += CG.flux_divergence_operator_from_advection(19e6u"m^3/s", [1, 3, 5, 6, 4, 2], v3D, nb)

    return T
end



"""
    load

Returns `grd` and `T` (in that order).
"""
function load()
    print("Creating 6-box model")
    grid = build_grid()
    T = build_T(grid)
    println(" ✔")
    @info """
            You are about to use the 3-box model of Archer et al. [2000].
            If you use it for research, please cite:

            - Archer, D. E., Eshel, G., Winguth, A., Broecker, W., Pierrehumbert, R., Tobis, M., and Jacob, R. (2000), Atmospheric pCO2 sensitivity to the biological pump in the ocean, Global Biogeochem. Cycles, 14 (4), 1219–1230, doi:10.1029/1999GB001216.

            You can find the corresponding BibTeX entry in the CITATION.bib file
            at the root of the AIBECS.jl package repository.
            (Look for the "Archer_etal_2000" key.)

            Note that although this model represents the 3-box model of Archer et al. [2000],
            it effectively requires 6 boxes to represent in the AIBECS, 
            so that particulate sinking flux divergence operators can be built.
            (See the comments at the start of the Archer_etal_2000.jl file for details.)
            """
    return grid, T
end

end

export Archer_etal_2000

