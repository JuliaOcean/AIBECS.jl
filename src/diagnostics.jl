


function stencil(T, grd::OceanGrid)
    # Figure out the stencil
    i, j, v = findnz(T) # linear indices of T
    # Cartesian indices on the full 3D space
    iwet = findall(vec(iswet(grd)))
    I = CartesianIndices(size(grd))[iwet[i]]
    J = CartesianIndices(size(grd))[iwet[j]]
    # Cartesian-index differences
    st = unique(J - I)
    # Taking the unique rows shows the stencil
    return unique(centereddiff(st, grd))
end
stencil(grd::OceanGrid, T) = stencil(T, grd)
export stencil

# Aux functions to center the stencil
function modcent(x, r)
    x = mod(x, r)
    return x < r/2 ? x : x - r
end
centereddiff(I, grd) = [CartesianIndex(modcent.(convert(Tuple,i), size(grd))) for i in I]





