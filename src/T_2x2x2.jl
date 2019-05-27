module T_2x2x2
#=
This circulation and grid was created by François Primeau (fprimeau@uci.edu)
It was taken from his Notebook at
https://github.com/fprimeau/BIOGEOCHEM_TEACHING/blob/master/Intro2TransportOperators.ipynb

The simple box model we consider is embeded in a 2×2×2 "shoebox".
It has 5 wet boxes and 3 dry boxes.
See image at
https://github.com/fprimeau/BIOGEOCHEM_TEACHING/blob/master/boxmodel.png
=#

using Unitful
using LinearAlgebra, SparseArrays
using ..GridTools

function build_wet3d()
    wet3d = trues(2, 2, 2)
    wet3d[[4,7,8]] .= false # land points
    return wet3d
end

# Macro to create grd in a simple `@Dict` call
# Thanks to Lyndon White (on slack)
macro Dict(vars...)
    kvs = Expr.(:call, :Pair, string.(vars), esc.(vars))
    Expr(:call, :Dict, kvs...)
end

function build_grd()
    wet3d = build_wet3d()
    # From F Primeau (slightly modified)
    earth_radius = 6367u"km"
    earth_area = 4π * earth_radius^2
    ocean_depth = 3700u"m"
    volume_ocean = earth_area * ocean_depth  # (including land points)
    h_top = 200u"m"        # thickness of top layer
    h_bottom = ocean_depth - h_top   # thickness of bottom layer
    # Build arrays for longitudes
    dxt = [180, 180]
    xt = [0, 180]
    XT3d = repeat(reshape(xt, (1,2,1)), outer=(2,1,2))
    DXT3d = repeat(reshape(dxt/360 * 2π * earth_radius, (1,2,1)), outer=(2,1,2)) # Careful: in meters here!
    # Build arrays for latitudes
    dyt = [90, 90]
    yt = [-45, 45]
    YT3d = repeat(reshape(yt, (2,1,1)), outer=(1,2,2))
    DYT3d = repeat(reshape(dyt/360 * 2π * earth_radius, (2,1,1)), outer=(1,2,2))
    # Build arrays for depths
    dzt = [h_top, h_bottom]
    zt = cumsum(dzt) - dzt/2
    zw = cumsum(dzt) - dzt
    ZT3d = repeat(reshape(zt, (1,1,2)), outer=(2,2,1))
    ZW3d = repeat(reshape(zw, (1,1,2)), outer=(2,2,1))
    DZT3d = repeat(reshape(dzt, (1,1,2)), outer=(2,2,1))
    return @Dict DXT3d DYT3d DZT3d ZW3d ZT3d dxt xt dyt yt dzt zt
end

function build_T(grd)
    wet3d = build_wet3d()
    iwet = LinearIndices(size(wet3d))[findall(wet3d)]
    v3d = grd["DXT3d"] .* grd["DYT3d"] .* grd["DZT3d"] .|> u"m^3"
    n = length(v3d)
    T = spzeros(n, n) * 1u"1/s"
    # Antarctic Circumpoloar Current 1 -> 3 -> 1 
    # TODO Add tools to GridTools to create those circulations and the mixings too
    ACC = 100e6u"m^3/s"
    isACC = sparse([1, 3], [3, 1], trues(2), n, n) ;
    # Meridional Overturning Circulation 1 -> 2 -> 6 -> 5 -> 1
    MOC = 15e6u"m^3/s"
    isMOC = sparse([1, 2, 6, 5], [2, 6, 5, 1], trues(4), n, n)
    # vertical mixing at "high northern latitudes"
    MIX = 10e6u"m^3/s"
    isMIX = sparse([2, 6], [6, 2], trues(2), n, n)
    # Fill in the divergences
    for orig in 1:n, dest in 1:n
        # Overturning part
        if isACC[orig, dest]
            T[orig, orig] += v3d[orig] \ ACC # add at origin
            T[dest, orig] -= v3d[dest] \ ACC # remove at destination
        end
        if isMOC[orig, dest]
            T[orig, orig] += v3d[orig] \ MOC # add at origin
            T[dest, orig] -= v3d[dest] \ MOC # remove at destination
        end
        if isMIX[orig, dest]
            T[orig, orig] += v3d[orig] \ MIX # add at origin
            T[dest, orig] -= v3d[dest] \ MIX # remove at destination
        end
    end
    return ustrip(T[iwet, iwet])
end

"""
    load

Returns wet3d, grd, and T (in that order).
"""
function load()
    print("Loading François Primeau's 2x2x2 model")
    wet3d = build_wet3d()
    grd = build_grd()
    T = build_T(grd)
    println(" ✅")
    return wet3d, grd, T
end

end

export T_2x2x2

