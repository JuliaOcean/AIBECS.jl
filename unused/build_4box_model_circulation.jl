const wet3d = trues(1, 2, 2)
using SparseArrays

# Macro to create grd in a simple `@Dict` call
# Thanks to Lyndon White (on slack)
macro Dict(vars...)
    kvs = Expr.(:call, :Pair, string.(vars), esc.(vars))
    Expr(:call, :Dict, kvs...)
end

function build_grd()
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
    DYT3d[1, 1, :] .= DYhigh
    DYT3d[1, 2, :] .= DYlow
    # Build DZ3d
    DZT3d = zeros(size(wet3d))
    DZT3d[1, 1, 1] = depth_high
    DZT3d[1, 2, 1] = depth_low
    DZT3d[1, 1, 2] = depth_ocean - depth_high
    DZT3d[1, 2, 2] = depth_ocean - depth_low
    # Build DX3d
    DXT3d = vtot / (depth_ocean * DYocean) * ones(size(wet3d))
    # ZT3d 
    ZT3d = cumsum(DZT3d, dims=3) - DZT3d / 2
    # ZW3d (depth at top of box)
    ZW3d = cumsum(DZT3d, dims=3) - DZT3d

    return @Dict DXT3d DYT3d DZT3d ZW3d ZT3d
end
const grd = build_grd()

function build_T(grd)
    # Mixing and overturning (F_HD, etc. see Figure 1 of Archer et al.)
    Mixing = zeros(4,4)
    Mixing[1, 2] = Mixing[2, 1] = 10e6 # between top boxes (1 and 2) (10 Sv = 10e6 cubic meters)
    Mixing[1, 3] = Mixing[3, 1] = 53e6 # between high-lat and deep (1 and 3)
    Mixing[2, 4] = Mixing[4, 2] = 1e6 # between low-lat and deep (2 and 4)
    Mixing[3, 4] = Mixing[4, 3] = 1e10 # (trick) fast mixing between deep boxes (3 and 4)
    overturning = 19e6
    Overturning = falses(4,4)
    Overturning[1, 3] = Overturning[3, 4] = Overturning[4, 2] = Overturning[2, 1] = true # Overturning circulation
    # Build T, as a divergence operator, i.e., 
    # T is positive where tracers are removed.
    T = spzeros(4,4)
    for orig in 1:4, dest in 1:4
        # Overturning part
        if Overturning[orig, dest] 
            T[orig, orig] += overturning # add at origin
            T[dest, orig] -= overturning # remove at destination
        end

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
const T = build_T(grd)

