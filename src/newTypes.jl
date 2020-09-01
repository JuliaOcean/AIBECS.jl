

"""
    ustrip(x::MetadataArray)

Strips unit from `x.parent` but stores it in `x.metadata` for safekeeping.
"""
Unitful.ustrip(x::MetadataArray) = MetadataArray(ustrip.(x.parent), (x.metadata..., unit=unit(eltype(x))))
"""
    upreferred(x::MetadataArray)

Converts `x.parent` to SI unit but keeps `x.metadata` for safekeeping.
"""
function Unitful.upreferred(x::MetadataArray)
    parent = upreferred.(x.parent)
    processing_str = "converted from $(unit(eltype(x.parent))) to preferred $(unit(eltype(parent)))"
    if haskey(x.metadata, :processing)
        metadata = x.metadata
        push!(metadata.processing, processing_str)
    else
        metadata = (x.metadata..., processing=[processing_str])
    end
    return MetadataArray(parent, metadata)
end
"""
    *(x::MetadataArray, q::Quantity)

    Preserve metadata (and history of unit conversions) when multiplying `x` by a quantity.
"""
function Base.:*(x::MetadataArray, q::Quantity)
    parent = x.parent * q
    processing_str = "× $q"
    if haskey(x.metadata, :processing)
        metadata = x.metadata
        push!(metadata.processing, processing_str)
    else
        metadata = (x.metadata..., processing=[processing_str])
    end
    return MetadataArray(parent, metadata)
end

"""
    onlykeep(x::MetadataArray, idx)

Returns `x[idx]` but also applies the index to the metadata that is originally of the same length as `x`.
"""
function onlykeep(x::MetadataArray, idx)
    xout = x.parent[idx]
    mout = (;((v isa Vector && length(v) == length(x)) ? (k,v[idx]) : (k,v) for (k,v) in pairs(x.metadata))...)
    return MetadataArray(xout, mout)
end
export onlykeep


"""
    AgedJacobianFactors

Type containing the Jacobian Factors and the age of the Jacobian.
This allows for the Shamanskii method to not update the Jacobian at each iterate.
"""
mutable struct AgedJacobianFactors{T}
    fac::T     # Jacobian factors — can be Real, Complex, Dual, or HyperDual
    age::Int   # age of the Jacobian
end

"""
    \\(Jf::JacobianFactors, y)

Dispatches backslash to work with all `JacobianFactors` subtypes.
"""
function Base.:\(Jf::AgedJacobianFactors, y)
    return Jf.fac \ y
end


