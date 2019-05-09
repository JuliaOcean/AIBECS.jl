# Reexport GridTools as they are useful outside too
@reexport module GridTools

using LinearAlgebra, SparseArrays
using WorldOceanAtlasTools
WOA = WorldOceanAtlasTools

#===================================
Off-diaongal matrices
===================================#

"""
    buildIbelow(wet3d, iwet)

Build the shifted-diagonal sparse matrix of the indices of below neighbours.

Ibelow[i,j] = 1 if the box represented by the linear index i
lies directly below the box represented by the linear index j.
"""
function buildIbelow(wet3d, iwet)
    nlat, nlon, ndepth = size(wet3d)
    n = nlon * nlat * (ndepth + 1)
    In = sparse(I, n, n)
    idx = zeros(Int64, nlat, nlon, ndepth + 1)
    idx[:] = 1:n
    idx .= idx[:, :, [2:ndepth + 1; 1]]      # downward shift
    return In[idx[:], :][iwet, iwet]
end
export buildIbelow

"""
    buildIabove(wet3d, iwet)

Build the shifted-diagonal sparse matrix of the indices of above neighbours.

Iabove[i,j] = 1 if the box represented by the linear index i
lies directly above the box represented by the linear index j.
"""
buildIabove(wet3d, iwet) = copy(transpose(buildIbelow(wet3d, iwet)))
export buildIabove



#===================================
Functions returning useful arrays,
vectors, and constants
===================================#

"""
    array_of_volumes(grd)

Returns the 3D array of volumes of grid boxes.
"""
array_of_volumes(grd) = grd["DXT3d"] .* grd["DYT3d"] .* grd["DZT3d"]
export array_of_volumes

"""
    indices_of_wet_boxes(wet3d)

Returns the vector of the indices of wet grid boxes.
"""
indices_of_wet_boxes(wet3d) = (LinearIndices(wet3d))[findall(!iszero, wet3d)]
export indices_of_wet_boxes

"""
    number_of_wet_boxes(wet3d)

Returns the number of wet grid boxes.
"""
number_of_wet_boxes(wet3d) = length(indices_of_wet_boxes(wet3d))
export number_of_wet_boxes

"""
    vector_of_volumes(wet3d, grd)

Returns the vector of volumes of wet boxes.
"""
vector_of_volumes(wet3d, grd) = array_of_volumes(grd)[indices_of_wet_boxes(wet3d)]
export vector_of_volumes

"""
    vector_of_depths(wet3d, grd)

Returns the vector of depths of the center of wet boxes.
"""
vector_of_depths(wet3d, grd) = grd["ZT3d"][indices_of_wet_boxes(wet3d)]
export vector_of_depths

"""
    vector_of_depths(wet3d, grd)

Returns the vector of depths of the top of wet boxes.
"""
vector_of_top_depths(wet3d, grd) = grd["ZW3d"][indices_of_wet_boxes(wet3d)]
export vector_of_top_depths

#===================================
Functions return weighted norms
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
Functions for the World Ocean Atlas
TODO
===================================#

"""
    observed_mean_and_variance(wet3d, grd, tracer;
                               product_year = 2018, # default WOA18
                               period = "Annual",   # default annual
                               resolution = "1°",   # default 1°×1° resolution
                               field = "an")        # default "objectively analyzed mean"

Returns the observed mean and variance of the tracer `tracer` fit to the `grd` grid.
AIBECS only uses the World Ocean Atlas for that at the moment.
Specifically, AIBECS uses the WorldOceanAtlasTools package.
However, this function will likely change since other observational datasets are available.
"""
function observed_mean_and_variance(wet3d, grd, tracer;
                                    product_year = 2018, # default WOA18
                                    period = "Annual",   # default annual
                                    resolution = "1°",   # default 1°×1° resolution
                                    field = "an")        # default "objectively analyzed mean"
    χ_3D, σ²_3D = WOA.fit_to_grid(grd, product_year, tracer, period, resolution, field)
    iwet = indices_of_wet_boxes(wet3d)
    return μ_3d[iwet], σ²_3d[iwet]
end
export observed_mean_and_variance

end

