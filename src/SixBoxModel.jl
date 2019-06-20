module SixBoxModel
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
using ..GridTools
using Unitful, UnitfulAstro # for units
using Reexport
@reexport using OceanGrids            # To store the grid

build_wet3D() = trues(2, 1, 3)

# Macro to create grid in a simple `@Dict` call
# Thanks to Lyndon White (on slack)
macro Dict(vars...)
    kvs = Expr.(:call, :Pair, string.(vars), esc.(vars))
    Expr(:call, :Dict, kvs...)
end

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
    return OceanGrid(elon, elat, edepth)
end


function build_T(grid)
    # Mixing and overturning (F_HD, etc. see Figure 1 of Archer et al.)
    Mixing = zeros(6, 6)
    Mixing[1, 2] = Mixing[2, 1] = 10e6 # between top boxes (1 and 2) (10 Sv = 10e6 cubic meters)
    Mixing[2, 5] = Mixing[5, 2] = 53e6 # between high-lat and deep (1 and 3)
    Mixing[2, 4] = Mixing[4, 2] = 1e6 # between low-lat and deep (2 and 4)
    # (trick) fast mixing between boxes of the same 3-box-model box
    idx_high = [1, 3]
    idx_low = [2]
    idx_deep = [4, 5, 6]
    for idx in [idx_high, idx_low, idx_deep]
        for i in idx, j in idx
            (i ≠ j) ? Mixing[i, j] = 1e10 : nothing
        end
    end
    # Overturning circulation
    overturning = 19e6
    Overturning = falses(6, 6)
    Overturning[1, 3] = Overturning[3, 5] = Overturning[5, 6] = Overturning[6, 4] = Overturning[4, 2] = Overturning[2, 1] = true 
    # Build T, as a divergence operator, i.e.,
    # T is positive where tracers are removed.
    T = spzeros(6,6)
    for orig in 1:6, dest in 1:6
        # Overturning part
        if Overturning[orig, dest]
            T[orig, orig] += overturning # add at origin
            T[dest, orig] -= overturning # remove at destination
        end
        # Mixing part
        if !iszero(Mixing[orig, dest])
            T[orig, orig] += Mixing[orig, dest] # add at origin
            T[dest, orig] -= Mixing[orig, dest] # remove at destination
        end
    end
    v3d = array_of_volumes(grid)
    v = vec(v3d) # Fine here because every point is wet :)
    V⁻¹ = sparse(Diagonal(v.^(-1)))
    T = V⁻¹ * T
    return T
end

"""
    load

Returns wet3D, grid, and T (in that order).
"""
function load()
    print("Loading 6-box model")
    wet3D = build_wet3D()
    grid = build_grid()
    T = build_T(grid)
    println(" ✅")
    return wet3D, grid, T
end

end

export SixBoxModel

