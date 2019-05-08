module SixBoxModelCirculation
#= 
The 6-box model is used because the code to bin WOA data
requires a regular grid.
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

using SparseArrays, SuiteSparse
using DataDeps, JLD2, FileIO

build_wet3d() = trues(2, 1, 3)

# Macro to create grd in a simple `@Dict` call
# Thanks to Lyndon White (on slack)
macro Dict(vars...)
    kvs = Expr.(:call, :Pair, string.(vars), esc.(vars))
    Expr(:call, :Dict, kvs...)
end

function build_grd()
    wet3d = build_wet3d()
    # From Archer et al. [2000]
    vtot = 1.292e18     # m³
    area_ocean = 349e12 # m²
    depth_low = 100.0  # m
    depth_high = 250.0 # m
    fraction_area_high = 0.15 #
    # Infer other values
    depth_ocean = vtot / area_ocean
    # Build DY3d
    earth_perimeter = 40e6 # 40'000 km
    DYocean = earth_perimeter / 2 # distance from south pole to north pole
    DYhigh = fraction_area_high * DYocean
    DYlow = (1 - fraction_area_high) * DYocean
    DYT3d = zeros(size(wet3d))
    DYT3d[1, :, :] .= DYhigh
    DYT3d[2, :, :] .= DYlow
    # Build DZ3d
    DZT3d = zeros(size(wet3d))
    DZT3d[:, :, 1] .= depth_low
    DZT3d[:, :, 2] .= depth_high - depth_low
    DZT3d[:, :, 3] .= depth_ocean - depth_high
    # Build DX3d
    DXT3d = vtot / (depth_ocean * DYocean) * ones(size(wet3d))
    # ZT3d
    ZT3d = cumsum(DZT3d, dims=3) - DZT3d / 2
    # ZW3d (depth at top of box)
    ZW3d = cumsum(DZT3d, dims=3) - DZT3d
    # lat, lon, depth
    dxt = [360.0]
    xt = dxt / 2
    dyt = [fraction_area_high * 180, (1 - fraction_area_high) * 180]
    yt = -90.0 .+ cumsum(dyt) .- dyt/2
    dzt = DZT3d[1, 1, :]
    zt = cumsum(dzt) .- dzt/2
    return @Dict DXT3d DYT3d DZT3d ZW3d ZT3d dxt xt dyt yt dzt zt
end

function build_T(grd)
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
    v3d = buildv3d(grd)
    v = vec(v3d) # Fine here because every point is wet :)
    V⁻¹ = d₀(v.^(-1))
    T = V⁻¹ * T
    return T
end

function load()
    println("loading 6-box-model circulation")
    wet3d = build_wet3d()
    grd = build_grd()
    T = build_T(grd)
    return wet3d, grd, T
end

end

