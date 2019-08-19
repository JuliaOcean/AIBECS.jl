# Reexport ObservationsTools as they are useful outside too
@reexport module ObservationsTools

using NearestNeighbors

"""
    interpolation_matrix(lat_obs, lon_obs, depth_obs, grd, wet3D)

Returns the matrix `M` such that `M * x` is the value of the tracer at obs locations.

The idea is to interpolate the modelled data onto the observed locations
rather than interpolating the observations onto the model grid.
For performant evaluation of the objective function,
it is maybe a good idea to create a constant matrix for the interpolation,
and to premultiply it in the expressions of the objective function and its derivative.
I.e., if

    f(x) = 1/2 δxᵀ W δx

where δx = x(model) - x(obs) with x(model) evaluated at the obs,

then x(model) = M * x for some matrix M (if nearest neighbors or linear interp.)
So

    f(x) = 1/2 (M x - obs)ᵀ W (M x - obs)

Importantly, 

    ∇f(x) = (M x - obs)ᵀ W M

So M is useful!
"""
function interpolation_matrix(lat_obs, lon_obs, depth_obs, grd, wet3D)
    iwet = findall(vec(wet3D))
    grdbox_centers = [ustrip.(grd.lat_3D[iwet]) ustrip.(grd.lon_3D[iwet]) ustrip.(grd.depth_3D[iwet])]
    grdbox_centers = permutedims(grdbox_centers, [2, 1])
    kdtree = KDTree(grdbox_centers)
    n_obs = length(lat_obs)
    I = 1:n_obs
    J = Int64[]
    for i in I
        push!(J, knn(kdtree, [lat_obs[i], lon_obs[i], depth_obs[i]], 1, true)[1][1])
    end
    V = trues(n_obs)
    n_wet = length(iwet)
    return sparse(I, J, V, n_obs, n_wet)
end

#=
Dummy data for testing

data = [ # lat     # lon    # depth    # value
          -9.86    145.0     200.0     1.0        ;
          -9.86    145.0     500.0     1.5        ;
          -9.86    145.0    1000.0     2.0        ;
          -9.86    145.0    2000.0     4.0        ;
          -9.86    165.0     200.0     2.0        ;
          -9.86    165.0     500.0     2.5        ;
          -9.86    165.0    1000.0     3.0        ;
          -9.86    165.0    2000.0     5.0        ;
         +15.32    145.0     200.0     1.0        ;
         +15.32    145.0     500.0     1.5        ;
         +15.32    145.0    1000.0     2.0        ;
         +15.32    145.0    2000.0     4.0        ;
         +15.32    165.0     200.0     2.0        ;
         +15.32    165.0     500.0     2.5        ;
         +15.32    165.0    1000.0     3.0        ;
       ]

lat_obs = data[:,1]
lon_obs = data[:,2]
depth_obs = data[:,3]
=#

end

