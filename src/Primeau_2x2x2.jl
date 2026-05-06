"""
    Primeau_2x2x2

Pedagogical 2×2×2 "shoebox" circulation by François Primeau, with 5 wet boxes
and 3 dry boxes (see his
[Intro2TransportOperators](https://github.com/fprimeau/BIOGEOCHEM_TEACHING/blob/master/Intro2TransportOperators.ipynb)
notebook and
[boxmodel.png](https://github.com/fprimeau/BIOGEOCHEM_TEACHING/blob/master/boxmodel.png)).
Use [`Primeau_2x2x2.load`](@ref) to obtain `(grd, T)`.
"""
module Primeau_2x2x2

using LinearAlgebra
using SparseArrays
using Unitful
using Reexport
@reexport using OceanGrids            # To store the grid
using ..CirculationGeneration
CG = CirculationGeneration


function build_wet3D()
    wet3D = trues(2, 2, 2)
    wet3D[[4, 7, 8]] .= false # land points
    return wet3D
end


function build_grid()
    elon = [0, 180, 360] * u"°"
    elat = [-90, 0, 90] * u"°"
    edepth = [0, 200, 3700] * u"m"
    return OceanGrid(elon, elat, edepth, build_wet3D())
end

function build_T(grid)
    # From Archer et al. [2000]
    v3D = array_of_volumes(grid)
    nb = length(v3D)

    # Antarctic Circumpoloar Current 1 -> 3 -> 1
    ACC = 100.0e6u"m^3/s"
    T = CG.T_advection(ACC, [1, 3, 1], v3D, nb)
    # Meridional Overturning Circulation 1 -> 2 -> 6 -> 5 -> 1
    MOC = 15.0e6u"m^3/s"
    T += CG.T_advection(MOC, [1, 2, 6, 5, 1], v3D, nb)
    # vertical mixing at "high northern latitudes" 2 <-> 6
    MIX = 10.0e6u"m^3/s"
    T += CG.T_diffusion(MIX, 2, 6, v3D, nb)

    # Only keep wet points
    iwet = findall(vec(grid.wet3D))

    return T[iwet, iwet]
end


"""
    load

Returns `grd` and `T` (in that order).
"""
function load()
    print("Creating François Primeau's 2x2x2 model")
    grd = build_grid()
    T = build_T(grd)
    println(" ✔")
    return grd, ustrip.(T)
end

end
