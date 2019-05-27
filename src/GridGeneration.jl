# Reexport GridTools as they are useful outside too
@reexport module GridGeneration

using LinearAlgebra, SparseArrays

#=============================================
Tools for making simple circulations
=============================================#

# Macro to create grd in a simple `@Dict` call
# Thanks to Lyndon White (on slack)
macro Dict(vars...)
    kvs = Expr.(:call, :Pair, string.(vars), esc.(vars))
    Expr(:call, :Dict, kvs...)
end

"""
    generate_grid(elon, elat, edepth)

Generates grid from rectilinear grid edges.
`elon`, `elat`, and `edepth` can have units.
If no unit is given, arguments are assumed 
to be degrees for latitudes and longitudes,
and meters for depths.
"""
generate_grid(elon::T, elat::T, edepth::T) where T = generate_grid(elon * u"°", elat * u"°", edepth * u"m")
function generate_grid(elon, elat, edepth, R=upreferred(6371.0u"km"))
    nx, ny, nz = length(elon) - 1, length(elat) - 1, length(edepth) - 1
    xt = 0.5 * (elon[1:end-1] + elon[2:end])
    yt = 0.5 * (elat[1:end-1] + elat[2:end])
    zt = 0.5 * (edepth[1:end-1] + edepth[2:end])
    dxt, dyt, dzt = diff(elon), diff(elat), diff(edepth)
    # depth objects
    zw = cumsum(dzt) - dzt
    ZT3d = repeat(reshape(zt, (1,1,nz)), outer=(ny,nx,1))
    ZW3d = repeat(reshape(zw, (1,1,nz)), outer=(ny,nx,1))
    DZT3d = repeat(reshape(dzt, (1,1,nz)), outer=(ny,nx,1))
    # lat objects
    YT3d = repeat(reshape(yt, (ny,1,1)), outer=(1,nx,nz))
    dyt_m = R * dyt ./ 360u"°"
    DYT3d = repeat(reshape(dyt_m, (ny,1,1)), outer=(1,nx,nz))
    # lon objects
    XT3d = repeat(reshape(xt, (1,nx,1)), outer=(ny,1,nz))
    # for DXT3d, I first calculate the area and then infer the equivalent distance in meters
    A2d = R^2 * abs.(sin.(elat[1:end-1]) - sin.(elat[2:end])) * ustrip.(dxt .|> u"rad")'
    A3d = repeat(A2d, outer=(1, 1, nz))
    DXT3d = A3d ./ DYT3d
    return @Dict DXT3d DYT3d DZT3d ZW3d ZT3d dxt xt dyt yt dzt zt
end

function flux_divergence_operator_from_advection(c, i, v3d, nb)
    unit(c[1] / v3d[1]) ≠ unit(u"1/s") && error("Unit problem")
    j = circshift(i,-1) # i = origin -> j = destination
    return sparse(i, i, ustrip(c ./ v3d[i]), nb, nb) - sparse(j, i, ustrip(c ./ v3d[j]), nb, nb)
end
function flux_divergence_operator_from_advection(c::Float64, i, v3d, nb)
    c = fill(c, length(i))
    flux_divergence_operator_from_advection(c, i, v3d, nb)
end    

end
