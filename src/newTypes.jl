
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



#= 
TODO
Work in progress to add some output types to use for plotting, etc.
=#

"""
    ModelOutput

Contains all the Model output.
To be sent to, e.g., plotting routines.
"""
struct ModelOutput
    name::String
    grid::OceanGrid
    tracers::Vector{Tracer}
    transports::Vector{Transport}
    fluxes::Vector{Flux}
    parameters::DataFrame
end

"""
    Tracer

Modelled tracer data and metadata
"""
struct Tracer
    name::String
    description::String
    values::Vector
    transport_name::String
end

"""
    Flux

Flux between two tracers
"""
struct Flux
    name::String
end

