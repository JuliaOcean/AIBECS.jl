


stencil(T, grd::OceanGrid) = unique(directions(T, grd))
stencil(grd::OceanGrid, T) = stencil(T, grd)
export stencil

# Aux functions to center the stencil
function modcent(x, r)
    x = mod(x, r)
    return x < r/2 ? x : x - r
end
function modcent(x::T, r::T) where {T<:CartesianIndex{3}}
    return CartesianIndex(modcent(x.I[1], r.I[1]), modcent(x.I[2], r.I[2]), modcent(x.I[3], r.I[3]))
end
function modcent!(x::Vector{T}, r::T) where {T<:CartesianIndex{3}}
    for (i, xᵢ) in enumerate(x)
        x[i] = modcent(xᵢ, r)
    end
    return x
end

function directions(T, grd)
    i, j, v = findnz(T) # linear indices of T
    return directions(i, j, grd)
end
function directions(i, j, grd)
    iwet = findall(vec(iswet(grd)))
    tmp = CartesianIndices(size(grd))
    @views JI = tmp[iwet[j]] - tmp[iwet[i]]
    modcent!(JI, CartesianIndex(size(grd)))
    return JI
end

# Labelling stencil directions
label(dir::CartesianIndex) = label(dir.I)
function label(dir::Tuple{T,T,T} where {T<:Int64})
    dy, dx, dz = dir
    return string(labelz(dz), labely(dy), labelx(dx))
end
labelz(dz::Int64) = dz == 0 ? "" : dz > 0 ? string("Below", labelz(dz - 1)) : string("Above", labelz(dz + 1))
labelx(dx::Int64) = dx == 0 ? "" : dx > 0 ? string("East", labelx(dx - 1)) : string("West", labelx(dx + 1))
labely(dy::Int64) = dy == 0 ? "" : dy > 0 ? string("North", labely(dy - 1)) : string("South", labely(dy + 1))

# Symmetric indices with values of T[i,j] and T[j,i]
function symmetric_IJVVᵀ(T)
    i, j, val = findnz(T)
    # symmetric indices
    ijsym = unique(zip([i;j], [j;i]))
    n = length(ijsym)
    I, J = Vector{Int64}(undef, n), Vector{Int64}(undef, n)
    val, valᵀ = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
    for (k, ij) in enumerate(ijsym)
        i, j = ij[1], ij[2]
        I[k], J[k] = i, j
        val[k], valᵀ[k] = T[i,j], T[j,i]
    end
    return I, J, val, valᵀ
end


function directional_transport(T, grd, dir)
    Isym, Jsym, Vsym, Vᵀsym = symmetric_IJVVᵀ(T)
    dirs = directions(Isym, Jsym, grd)
    n = size(T,1)
    v = volumevec(grd)
    return _directional_transport(dir, dirs, n, Isym, Jsym, Vsym, Vᵀsym, v)
end

function _directional_transport(dir, dirs, n, Isym, Jsym, Vsym, Vᵀsym, v)
    all(Tuple(dir) .== 0) && error("No self directional transport")
    idir = findall(isequal(dir), dirs)
    @views i, j, val, valᵀ = Isym[idir], Jsym[idir], Vsym[idir], Vᵀsym[idir]
    return sparse([i; i], [j; i], [val; -v[j] ./ v[i] .* valᵀ], n, n)
end
function directional_transports(T, grd)
    st = stencil(grd, T)
    Isym, Jsym, Vsym, Vᵀsym = symmetric_IJVVᵀ(T)
    dirs = directions(Isym, Jsym, grd)
    n = size(T,1)
    v = volumevec(grd)
    return Dict((label(dir), _directional_transport(dir, dirs, n, Isym, Jsym, Vsym, Vᵀsym, v)) for dir in st if !isempty(label(dir)))
end
export directional_transport, directional_transports







