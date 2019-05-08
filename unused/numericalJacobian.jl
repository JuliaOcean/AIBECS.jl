# Numerical Jacobian for a sparsity pattern thta only involves neighbors
using SparseArrays
"""
    sparse(IJ::Vector{<:CartesianIndex}, v, m, n)

Overloads SparseArrays' `sparse` function to work directly with `CartesianIndex` type of index.
"""
function SparseArrays.sparse(IJ::Vector{<:CartesianIndex}, v, m, n)
    IJ′ = reinterpret(Int, reshape(IJ, 1, :))
    return sparse(view(IJ′, 1, :), view(IJ′, 2, :), v, m, n)
end

using DualNumbers
"""
    Compressor

Type, containing all the required objects for fast comoression and decompression.
"""
struct Compressor
    i
    j
    icmp
    jcmp
    ijcmp
    idcmp
    jdcmp
    nwet
    ncol
    ε_colors
    nz
    col_idx
end

"""
    compress(Jac, C::Compressor)

Compresses the Jacobian into a smaller matrix.
"""
function compress(Jac, C::Compressor)
    out = zeros(eltype(Jac), C.nwet, C.ncol)
    out[C.ijcmp] .= Jac.nzval
    return out
end

"""
    decompress(Jac, C::Compressor)

Decompresses the Jacobian back into its full size.
"""
decompress(smallJac, C::Compressor) = sparse(C.i, C.j, smallJac[C.ijcmp], C.nwet, C.nwet)

"""
    numJac(f, x, C)

Evaluates the Jacobian of `f` at `x` numerically using the `DualNumbers` package.

The evaluation is efficient by limiting the number of evaluations of `f`.
"""
function numJac(f, x, C::Compressor)
    xpert = repeat(x, 1, C.ncol) + C.ε_colors
    smallJac = epsilon.(f(xpert))
    J = decompress(smallJac, C)
    dropzeros!(J)
    return J
end

"""
    adjacency_matrix_indices(wet3d, dist)

Builds the adjacency matrix indices for the full 3d grid (wet and dry).

`dist` is the infitiny-norm distance of adjacent boxes.
Color indices thus vary from `1` to `(1+2dist)^2` boxes are adjacent for each box.
(`dist = 2` required for OCIM1.1 and POP2 and means 125 adjacent boxes.)
"""
function adjacency_matrix_indices(wet3d, dist)
    nlat, nlon, ndepth = size(wet3d)
    L = LinearIndices(wet3d)
    R = CartesianIndices(size(wet3d)) # Range of all cartesian indices
    i1 = first(R) # convenient to create R1 for accessing neighbors in the line below
    RAdj = CartesianIndices(UnitRange.(Tuple(-dist * i1), Tuple(dist * i1))) # Range of cartesian indices of neighbors 3 directions
    L1 = LinearIndices(size(RAdj))
    nR, nRAdj = length(R), length(RAdj)
    # Pre-allocate sub indices
    ljnz = zeros(Int64, nR, nRAdj)
    linz = zeros(Int64, nR)
    for i in R # Loop through all boxes
        li = L[i]
        linz[li] = li # Allocate index of each box
        for k in L1 # Loop through all adjacent boxes
            lj = L[valid_index(i + RAdj[k], nlat, nlon, ndepth)]
            ljnz[L[i], L1[k]] = lj # Allocate index of each adjacent box
        end
    end
    return linz, ljnz
end

"""
    adjacency_matrix(wet3d)

Builds the adjacency matrix for the full 3d grid (wet and dry).

`dist` is the infitiny-norm distance of adjacent boxes.
Color indices thus vary from `1` to `(1+2dist)^2` boxes are adjacent for each box.
(`dist = 2` required for OCIM1.1 and POP2 and means 125 adjacent boxes.)
"""
function adjacency_matrix(wet3d, dist)
    linz, ljnz = adjacency_matrix_indices(wet3d, dist)
    linz = repeat(linz, size(ljnz, 2))
    ljnz = vec(ljnz)
    AM_3d = sparse(linz, ljnz, true)
    return AM_3d
end

"""
    valid_index(j, nlat, nlon, ndepth)

Return a valid cartesian index.

If out of bounds in latitude or depth, return the bound.
If out of bounds in longitude, return the modulo(nlon) value.
(Because longitude is periodic.)
"""
function valid_index(j, nlat, nlon, ndepth)
    y = max(1, min(nlat, j[1]))
    x = mod1(j[2], nlon)
    z = max(1, min(ndepth, j[3]))
    return CartesianIndex(y, x, z)
end

"""
    build_compressor(wet3d, iwet, nwet, dist)

Builds a Compressor object.

A Compressor object contains an adjacency matrix and all the indices required
for compressing and decompressing the numerical Jacobian.
"""
function build_compressor(wet3d, iwet, nwet, dist)
    # Decompression indices
    # We use the indices of the adjacency matrix
    AM_3d = adjacency_matrix(wet3d, dist) # Full wet and dry adjacency adjacency_matrix
    AM = AM_3d[iwet, iwet] # Adjacency matrix (wet points only)
    i, j, _ = findnz(AM) # Sub indices

    # Compression indices
    # We use the color indices
    # The compressed index `icmp` is the same as the decompressed `i`
    # The compressed index `jcmp` is the color index
    icmp = i
    ncol, ε_colors, col_idx = color_groups(wet3d, iwet, nwet, dist)
    jcmp = col_idx[j]

    # matrix for preallocating smallJac (for storage in Compressor)
    #smallJac₀ = sparse(icmp, jcmp, 1.0, nwet, ncol)
    ijcmp = CartesianIndex.(icmp, jcmp)

    # number of non zeros for storage and use in numJac
    nz = length(i)

    # Permuted decompression indices (faster!)
    v = 1:nz # to track the index permutation that happense during compression
    AMv = sparse(icmp, jcmp, v, nwet, ncol)
    _, _, v2 = findnz(AMv) # v2 is the verctor of the index permutation
    idcmp, jdcmp = i[v2], j[v2] # This allows direct use of smallJac.nzval
    # which is much faster than accessing those values via (icmp, jcmp)

    return Compressor(i, j, icmp, jcmp, ijcmp, idcmp, jdcmp, nwet, ncol, ε_colors, nz, col_idx)
end

"""
    color_groups(wet3d, iwet, nwet, dist)

Builds color groups.
Outputs:
ncol, ε_colors, col_idx
`ncol` is the number of colors (should be (1+2dist)^2)
`ε_colors`` is the matrix of 1.0's used to perturb `x` in all directions of each color.
`col_idx` is the vector of the colors of each box.

Inputs:
`dist` is the infitiny-norm distance of adjacent boxes.
(`dist = 2` required for OCIM1.1 and POP2 and means 125 adjacent boxes.)
"""
function color_groups(wet3d, iwet, nwet, dist)
    nlon = size(wet3d, 2)
    if (nlon % (2dist + 1) ≠ 0) & (nlon ≠ 1)
        error("nlon must be a multiple of 2*dist+1!!")
    end
    col_idx_3d = zeros(Int64, size(wet3d))
    color = 0
    n = 1 + 2dist # dist is the maximum distance of adjacent boxes
    for idepth in 1:n, ilon in 1:n, ilat in 1:n
        color += 1
        # @show (ilat, ilon, idepth)
        col_idx_3d[ilat:n:end, ilon:n:end, idepth:n:end] .= color
    end
    nlon == 1 ? col_idx_3d[:] .= collect(1:nwet) : nothing # hack for 4-box-model
    col_idx = col_idx_3d[iwet]

    ncol = n^3 # number of colors
    ε_colors = zeros(nwet, ncol) # for perturbation in numJac
    for color = 1:ncol
        indices_of_color = findall(col_idx .== color)
        ε_colors[indices_of_color, color] .= 1.0
    end
    ε_colors = ε * ε_colors # make it dual-unit valued

    return ncol, ε_colors, col_idx
end

function build_big_compressor(wet3d, iwet, nwet, dist, nTracers)
    # Decompression indices
    # We use the indices of the adjacency matrix
    AM_3d = adjacency_matrix(wet3d, dist) # Full wet and dry adjacency adjacency_matrix
    AM = AM_3d[iwet, iwet] # Adjacency matrix (wet points only)
    i, j, _ = findnz(AM) # Sub indices

    # Compression indices
    # We use the color indices
    # The compressed index `icmp` is the same as the decompressed `i`
    # The compressed index `jcmp` is the color index
    icmp = i
    ncol, ε_colors, col_idx = color_groups(wet3d, iwet, nwet, dist)
    jcmp = col_idx[j]

    # number of non zeros for storage and use in numJac
    nz = length(i)

    # make everything size nTracers x nTracers
    bigAM = kron(sparse(trues(nTracers, nTracers)), AM)
    bigi, bigj, _ = findnz(bigAM) # Sub indices

    big_col_idx = zeros(eltype(col_idx), nTracers * nwet)
    for t in 0:nTracers-1
        big_col_idx[(1:nwet) .+ t * nwet] .= col_idx .+ t * ncol
    end

    bigicmp = bigi
    bigjcmp = big_col_idx[bigj]
    bigijcmp = CartesianIndex.(bigicmp, bigjcmp)

    bigε_colors = kron(Matrix(I, nTracers, nTracers), ε_colors)

    # Permuted decompression indices (faster!)
    bigncol = ncol * nTracers
    bignwet = nwet * nTracers
    bigv = 1:(nz * nTracers^2) # to track the index permutation that happense during compression
    bigAMv = sparse(bigicmp, bigjcmp, bigv, bignwet, bigncol)
    bigicmp2, bigjcmp2, bigv2 = findnz(bigAMv) # v2 is the verctor of the index permutation
    bigidcmp, bigjdcmp = bigi[bigv2], bigj[bigv2] # This allows direct use of smallJac.nzval
    # which is much faster than accessing those values via (icmp, jcmp)


    return Compressor(bigi, bigj, bigicmp, bigjcmp, bigijcmp, bigidcmp, bigjdcmp, bignwet, bigncol, bigε_colors, nz, big_col_idx)
end