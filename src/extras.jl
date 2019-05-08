#===================================
Constants
===================================#

"""
    constants(wet3d, grd)

Returns `nb`, `DIV`, `Iabove`, `v`, `z`, `ztop`.
These are:
- `nb` — the number of wet grid boxes
- `DIV` — the divergence operator
- `Iabove` — the matrix shifting indices to grid box above
- `v` — the vector of volumes of each grid box
- `z` — the vector of depths at the center of each grid box
- `ztop` — the vector of depths at the top of each grid box
These constants are packaged together in this function for compactness.
Among other things, they can be used to create:
- the multi-tracer model formulation
- the transport matrices for sinking particles (particle flux divergences)
"""
function constants(wet3d, grd)
    iwet = (LinearIndices(wet3d))[findall(!iszero, wet3d)]
    nb = length(iwet)
    DIV = buildDIV(wet3d, iwet, grd)
    Iabove = buildIabove(wet3d, iwet)
    v3d = buildv3d(grd)
    v = v3d[iwet]
    z = grd["ZT3d"][iwet]
    ztop = grd["ZW3d"][iwet]
    return nb, DIV, Iabove, v, z, ztop
end
export constants

#===================================
Functions: norms
===================================#

"""
    vnorm²(x, v)

Returns the square of the volume-weighted norm of tracer `x` using volumes `v`.
"""
vnorm²(x, v) = transpose(x) * Diagonal(v) * x

"""
    vnorm(x, v)

Returns the volume-weighted norm of tracer `x` using volumes `v`.
"""
vnorm(x, v) = sqrt(vnorm²(x, v))

"""
    vmean(x, v)

Returns the volume-weighted mean of tracer `x` using volumes `v`.
"""
vmean(x, v) = transpose(v) * x / sum(v)

"""
    vmean(x, v, i)

Returns the volume-weighted mean of tracer `x` using volumes `v`, but only over indices `i`.
This is useful to get the observed mean.
(Because there are some grid boxes without observations.)
"""
vmean(x, v, i) = transpose(v[i]) * x[i] / sum(v[i])

#===================================
Functions: World Ocean Atlas
===================================#

function build_μ_and_σ²_from_WOA(wet3d, grd, vv)
    iwet = (LinearIndices(wet3d))[findall(!iszero, wet3d)]
    μ_3d, σ²_3d = WorldOceanAtlasTools.WOA13_mean_and_variance_per_grid_box(grd, vv, "Annual", "1°")
    return μ_3d[iwet], σ²_3d[iwet]
end
