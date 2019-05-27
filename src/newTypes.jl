
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


