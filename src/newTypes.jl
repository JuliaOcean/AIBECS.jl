

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
    processing_str = "converted to preferred $(unit(eltype(parent)))"
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
    processing_str = "converted to $(unit(eltype(parent))) with $q"
    if haskey(x.metadata, :processing)
        metadata = x.metadata
        push!(metadata.processing, processing_str)
    else
        metadata = (x.metadata..., processing=[processing_str])
    end
    return MetadataArray(parent, metadata)
end

"""
    AgedJacobianFactors

Type containing the Jacobian Factors and the age of the Jacobian.
This allows for the Shamanskii method to not update the Jacobian at each iterate.
"""
mutable struct AgedJacobianFactors
    fac        # Jacobian factors — can be Real, Complex, Dual, or HyperDual
    age::Int # age of the Jacobian
    facage::Int # age of the factors of the Jacobian
end

"""
    \\(Jf::JacobianFactors, y)

Dispatches backslash to work with all `JacobianFactors` subtypes.
"""
function Base.:\(Jf::AgedJacobianFactors, y)
    return Jf.fac \ y
end


