
abstract type SubSparsityPattern end
abstract type TransportPattern <: SubSparsityPattern end
abstract type TransferPattern <: SubSparsityPattern end
struct CirculationSparsity <: TransportPattern end
struct ParticleSparsity <: TransportPattern end
struct DiagonalSparsity <: TransferPattern end
struct ZeroSparsity <: TransferPattern end

export CirculationSparsity, ParticleSparsity, DiagonalSparsity, ZeroSparsity


Base.show(io::IO, ::MIME"text/plain", ::DiagonalSparsity) = print(io, "⠑⢄")
Base.show(io::IO, ::MIME"text/plain", ::CirculationSparsity) = print("⢿⣷")
Base.show(io::IO, ::MIME"text/plain", ::ParticleSparsity) = print("⠳⣄")
Base.show(io::IO, ::MIME"text/plain", ::ZeroSparsity) = print("  ")
function Base.show(io::IO, m::MIME"text/plain", M::Array{<:SubSparsityPattern,2})
    for r in eachrow(M)
        for x in r
            show(io, m, x)
        end
        print("\n")
    end
end


#function sparseones(s::CirculationSparsity, grd)
#    st = s.st
#    LIA = LinearIndices(size(grd))
#    CIA = CartesianIndices(size(grd))
#    i = collect(LIA)
#    J = sum(sparse(i, LIA[CIA .+ dir], v, n, n) for dir in stencil(T, grd))
#end
#
#
#function grdindices(ijk::CartesianIndex{3}, grd)
#    i, j, k = Tuple(ijk)
#    CartesianIndex()
#end

import SparseArrays.sparse
function sparse(::DiagonalSparsity, nb, T)
    sparse(Diagonal(trues(nb)))
end
function sparse(::CirculationSparsity, nb, T)
    i, j = findnz(T)
    sparse(i, j, trues(length(i)), nb, nb)
end
function sparse(::ParticleSparsity, nb, T)
    tmp = buildIabove(grd) + I # <- this does not work because no grd available here!!!
    i, j = findnz(tmp)
    sparse(i, j, trues(length(i)), nb, nb)
end
function sparse(::ZeroSparsity, nb, T)
    sparse([], [], trues(0), nb, nb)
end
function sparse(M::Array{<:SubSparsityPattern,2}, nb, T)
    m, n = size(M)
    m ≠ n && error("Incorrect size")
    reduce(vcat, reduce(hcat, sparse(M[i,j], nb, T) for j in 1:n) for i in 1:m)
end
export sparse


"""
    findshiftedblockindices(jac, nb, nt, i, j)

Returns the shifted `(I, J)` indices of `jac` of the `(i, j)`-th block
"""
function findshiftedblockindices(jac, nb, nt, i, j)
    m, n = size(jac)
    ((m ≠ n) || (nb * nt ≠ n)) && error("Incorrect sizes")
    I, J = findnz(jac)
    I .-= (i-1)*nb
    J .-= (j-1)*nb
    I, J
end
"""
    findblockindices(jac, nb, nt, i, j)

Returns the indices of `jac.nzval` of the `(i, j)`-th block of `jac`.
To be used for in-place sparse Jacobian updates.
"""
function findblockindices(jac, nb, nt, i, j)
    I, J = findshiftedblockindices(jac, nb, nt, i, j)
    findall(@. (0 < I ≤ nb) & (0 < J ≤ nb))
end
"""
    findblockdiagonalindices(jac, nb, nt, i, j)

Returns the indices of `jac.nzval` of the diagonal of the `(i, j)`-th block of `jac`.
To be used for in-place sparse Jacobian updates.
(This is different from `findblockindices` and must be used to fill in the diagonal terms
of diagonal blocks, which are shared between the Jacobians of the `G`s and the `T`s.)
"""
function findblockdiagonalindices(jac, nb, nt, i, j)
    I, J = findshiftedblockindices(jac, nb, nt, i, j)
    findall(@. (0 < I == J ≤ nb))
end

