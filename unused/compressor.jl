function compress(X, i_c, j_c, ij_d)
# Compresses a sparse, square matrix X into a rectangular matrix where each column is a full layer
    X_comp = sparse(i_c, j_c, X[ij_d])
    return full(X_comp)
end

function decompress(X_comp, i_d, j_d, ij_c)
# Decompresses a rectangular compressed matrix X into
# a square sparse matrix.
    X = sparse(i_d, j_d, X_comp[ij_c])
    return X
end

# TODO Make sure I can use these tricks to compute a quick analytical Jacobian
# for use in ∂xregen!
function waterColumnStructure(wet3d)
    Iabove = buildIabove(wet3d)
  # loop to make all powers of A
    B = spones(Iabove + I)
    foo = true
    while foo
        last_B = copy(B)
        B = spones(B^2)
        foo = B != last_B
    end
    return B
end

function buildCompressionIndices(wet3d)
  # TODO describe what this does
    layer_idx_3d = zeros(wet3d, Int64)
    ndepth = size(wet3d, 3)
    iwet = find(wet3d)
    for idepth = 1:ndepth
        layer_idx_3d[:,:,idepth] .= idepth
    end
    layer_idx = layer_idx_3d[iwet]

    B = waterColumnStructure(wet3d)
    ij_d = find(B)
    i_d, j_d = ind2sub(size(B), ij_d)

    i_c = copy(i_d)
    j_c = layer_idx[j_d] # This is because I know that this is the right scheme
    B_c = compress(B, i_c, j_c, ij_d)
    ij_c = sub2ind(size(B_c), i_c, j_c)
    return i_c, j_c, ij_c, i_d, j_d, ij_d
end

export buildCompressionIndices, compress, decompress
