# Reexport GridTools as they are useful outside too
@reexport module GridTools

using LinearAlgebra, SparseArrays

#===================================
Off-diaongal matrices
===================================#

"""
    buildIbelow(wet3D, iwet)

Builds the shifted-diagonal sparse matrix of the indices of below neighbours.

Ibelow[i,j] = 1 if the box represented by the linear index i
lies directly below the box represented by the linear index j.
"""
function buildIbelow(wet3D, iwet)
    nlat, nlon, ndepth = size(wet3D)
    n = nlon * nlat * (ndepth + 1)
    In = sparse(I, n, n)
    idx = zeros(Int, nlat, nlon, ndepth + 1)
    idx[:] = 1:n
    idx .= idx[:, :, [2:ndepth + 1; 1]]      # downward shift
    return In[idx[:], :][iwet, iwet]
end
"""
    buildIbelow(grid)

Builds the shifted-diagonal sparse matrix of the indices of below neighbours for `grid`.
See `buildIbelow(wet3D, iwet)`.
"""
buildIbelow(grid) = buildIbelow(grid.wet3D, indices_of_wet_boxes(grid))
export buildIbelow

"""
    buildIabove(wet3D, iwet)

Builds the shifted-diagonal sparse matrix of the indices of above neighbours.

Iabove[i,j] = 1 if the box represented by the linear index i
lies directly above the box represented by the linear index j.
"""
buildIabove(wet3D, iwet) = copy(transpose(buildIbelow(wet3D, iwet)))
"""
    buildIabove(grid)

Builds the shifted-diagonal sparse matrix of the indices of above neighbours for `grid`.
See `buildIabove(wet3D, iwet)`.
"""
buildIabove(grid) = copy(transpose(buildIbelow(grid)))
export buildIabove



#===================================
Functions returning useful arrays,
vectors, and constants
===================================#

"""
    array_of_volumes(grid)

Returns the 3D array of volumes of grid boxes.
"""
array_of_volumes(grid) = grid.volume_3D
export array_of_volumes


"""
    indices_of_wet_boxes(wet3D::BitArray)

Returns the vector of the indices of wet grid boxes.
"""
indices_of_wet_boxes(wet3D::BitArray) = findall(vec(wet3D))
"""
    indices_of_wet_boxes(grid)

Returns the vector of the indices of wet grid boxes.
"""
indices_of_wet_boxes(grid) = indices_of_wet_boxes(grid.wet3D)
export indices_of_wet_boxes

"""
    number_of_wet_boxes(wet3D::BitArray)

Returns the number of wet grid boxes.
"""
number_of_wet_boxes(wet3D::BitArray) = length(indices_of_wet_boxes(wet3D))
"""
    number_of_wet_boxes(grid)

Returns the number of wet grid boxes.
"""
number_of_wet_boxes(grid) = number_of_wet_boxes(grid.wet3D)
export number_of_wet_boxes

"""
    vector_of_volumes(grid)

Returns the vector of volumes of wet boxes.
"""
vector_of_volumes(grid) = array_of_volumes(grid)[indices_of_wet_boxes(grid.wet3D)]
export vector_of_volumes

"""
    vector_of_depths(grid)

Returns the vector of depths of the center of wet boxes.
"""
vector_of_depths(grid) = grid.depth_3D[indices_of_wet_boxes(grid.wet3D)]
export vector_of_depths

"""
    vector_of_depths(grid)

Returns the vector of depths of the top of wet boxes.
"""
vector_of_top_depths(grid) = grid.depth_top_3D[indices_of_wet_boxes(grid.wet3D)]
export vector_of_top_depths

#===================================
Functions return weighted norms
===================================#

"""
    weighted_norm²(x, w)

Returns the square of the weighted norm of `x` using weights `w`.
"""
weighted_norm²(x, w) = transpose(x) * Diagonal(w) * x

"""
    weighted_norm(x, v)

Returns the weighted norm of `x` using weights `w`.
"""
weighted_norm(x, w) = sqrt(weighted_norm²(x, w))

"""
    weighted_mean(x, w)

Returns the weighted mean of `x` using weights `w`.
"""
weighted_mean(x, w) = transpose(w) * x / sum(w)

"""
    vmean(x, w, I)

Returns the weighted mean of `x` using weights `w`, but only over indices `I`.
This is useful to get the observed mean.
(Because there are some grid boxes without observations.)
"""
weighted_mean(x, w, I) = transpose(w[I]) * x[I] / sum(w[I])
export weighted_norm²

#===================================
1D Vector <-> 3D array conversions
===================================#

"""
    rearrange_into_3Darray(x, wet3D::BitArray)

Returns a 3D array of `x` rearranged to the `true` entries of `wet3D`.
Entries where `wet3D` is `false` are filled with `NaN`s.
"""
function rearrange_into_3Darray(x, wet3D::BitArray)
    iwet = indices_of_wet_boxes(wet3D)
    x3d = Array{Union{Missing,Float64},3}(undef, size(wet3D))
    x3d[iwet] .= x
    return x3d
end
"""
    rearrange_into_3Darray(x, grid)

Returns a 3D array of `x` rearranged to the wet boxes of the grid.
"""
rearrange_into_3Darray(x, grid) = rearrange_into_3Darray(x, grid.wet3D)

export rearrange_into_3Darray

#=============================================
unpacking of multi-tracers
=============================================#

state_to_tracers(x, nb, nt) = ntuple(i -> state_to_tracer(x, nb, nt, i), nt)
state_to_tracer(x, nb, nt, i) = x[tracer_indices(nb, nt, i)]
function state_to_tracers(x, grd) 
    nb = number_of_wet_boxes(grd)
    nt = Int(round(length(x) / nb))
    return state_to_tracers(x, nb, nt)
end
tracer_indices(nb, nt, i) = (i-1)*nb+1 : i*nb
tracers_to_state(xs) = reduce(vcat, xs)
export state_to_tracers, state_to_tracer, tracers_to_state, tracer_indices
# Alias for better name
unpack_tracers = state_to_tracers
export unpack_tracers

end

