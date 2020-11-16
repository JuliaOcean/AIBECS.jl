import Flatten: flattenable
@metadata initial_value nothing
import FieldMetadata: @units, units
import FieldMetadata: @limits, limits
import FieldMetadata: @prior, prior
import FieldMetadata: @description, description
@metadata bounds nothing
import FieldMetadata: @logscaled, logscaled
@metadata reference nothing


"""
    AbstractParameters{T} <: AbstractVector{T}

An abstract type for AIBECS model parameters.

Parameters in AIBECS use the following convenience packages:

- Parameters
- FieldMetadata
- FieldDefaults
- Flatten
- Unitful
- DataFrames
- Distributions

These aim to allow for some nice features, which include

- nice syntax for unpacking parameters in functions via the `@unpack` macro (fron UnPack.jl)
- additional metadata on parameters
- easy conversion to and from vectors
- use of units and automatic conversions if necessary
- pretty table-format displays
- loading and saving to and from CSV files
- prior estimates for bayesian inference and optimization

See the list of examples to get an idea of how to generate parameters for your model.

# Examples

Generate a simple parameter type via

```jldoctest
julia> struct SimpleParams{T} <: AbstractParameters{T}
           α::T
           β::T
           γ::T
       end
SimpleParams
```

To create an instance of the `SimpleParams(Float64)` type, you can do

```jldoctest
julia> p = SimpleParams(1.0, 2.0, 3.0)
SimpleParams{Float64}
│ Row │ Symbol │ Value   │
│     │ Symbol │ Float64 │
├─────┼────────┼─────────┤
│ 1   │ α      │ 1.0     │
│ 2   │ β      │ 2.0     │
│ 3   │ γ      │ 3.0     │
```

One of the core features from Parameters is unpacking in functions, e.g.,

```jldoctest
julia> function simplef(p)
           @unpack α, γ = p
           return α + γ
       end
simplef (generic function with 1 method)

julia> simplef(p) # 1.0 + 3.0
4.0
```

More complex examples are permitted by adding metadata (thanks to FieldMetadata.jl).
You can add units

```jldoctest
julia> @units struct UnitParams{T} <: AbstractParameters{T}
           α::T | u"km"
           β::T | u"hr"
           γ::T | u"m/s"
       end ;

julia> p = UnitParams(1.0, 2.0, 3.0)
UnitParams{Float64}
│ Row │ Symbol │ Value   │ Unit     │
│     │ Symbol │ Float64 │ Unitful… │
├─────┼────────┼─────────┼──────────┤
│ 1   │ α      │ 1.0     │ km       │
│ 2   │ β      │ 2.0     │ hr       │
│ 3   │ γ      │ 3.0     │ m s^-1   │
```

Note that when adding units to your parameters, they will be converted to SI
when unpacked, as in, e.g.,

```jldoctest
julia> function speed(p)
           @unpack α, β, γ = p
           return α / β + γ
       end
speed (generic function with 1 method)

julia> speed(p) # (1.0 km / 2.0 hr + 3 m/s) in m/s
3.138888888888889
```

Another example for optimizable/flattenable parameters

```jldoctest
julia> @initial_value @units @flattenable struct OptParams{T} <: AbstractParameters{T}
           α::T | 3.6 | u"km"  | true
           β::T | 1.0 | u"hr"  | false
           γ::T | 1.0 | u"m/s" | true
       end ;

julia> p = OptParams(initial_value(OptParams)...)
OptParams{Float64}
│ Row │ Symbol │ Value   │ Initial value │ Unit     │ Optimizable │
│     │ Symbol │ Float64 │ Float64       │ Unitful… │ Bool        │
├─────┼────────┼─────────┼───────────────┼──────────┼─────────────┤
│ 1   │ α      │ 3.6     │ 3.6           │ km       │ 1           │
│ 2   │ β      │ 1.0     │ 1.0           │ hr       │ 0           │
│ 3   │ γ      │ 1.0     │ 1.0           │ m s^-1   │ 1           │
```

Thanks to the FieldMetaData interface, you can chain the following preloaded metadata:

- initial_value
- units (from Unitful.jl)
- prior (from Distributions.jl)
- description (`String`)
- bounds (2-element `Tuple`)
- logscaled (`Bool`)
- flattenable (to convert to vectors of optimizable parameters only)
- reference (`String`)

Here is an example of parameter with all the possible metadata available in AIBECS:

```jldoctest
julia> @initial_value @units @prior @description @bounds @logscaled @flattenable @reference struct FullParams{T} <: AbstractParameters{T}
           α::T | 1.0 | u"km"  | Normal(0,1)    | "The distance"   | (-Inf, Inf) | false | false | "Jean et al., 2042"
           β::T | 2.0 | u"hr"  | LogNormal(0,1) | "The time"       | (   0, Inf) | true  | true  | "Claude et al. 1983"
           γ::T | 3.0 | u"mol" | Normal(1,2)    | "The # of moles" | (  -1,   1) | false | true  | "Dusse et al. 2000"
       end ;

julia> FullParams(4.0, 5.0, 6.0)
FullParams{Float64}
│ Row │ Symbol │ Value   │ Initial value │ Unit     │ Prior                            │ Description    │ Bounds      │ Logscaled │ Optimizable │ Reference          │
│     │ Symbol │ Float64 │ Float64       │ Unitful… │ Distribu…                        │ String         │ Tuple…      │ Bool      │ Bool        │ String             │
├─────┼────────┼─────────┼───────────────┼──────────┼──────────────────────────────────┼────────────────┼─────────────┼───────────┼─────────────┼────────────────────┤
│ 1   │ α      │ 4.0     │ 1.0           │ km       │ Normal{Float64}(μ=0.0, σ=1.0)    │ The distance   │ (-Inf, Inf) │ 0         │ 0           │ Jean et al., 2042  │
│ 2   │ β      │ 5.0     │ 2.0           │ hr       │ LogNormal{Float64}(μ=0.0, σ=1.0) │ The time       │ (0, Inf)    │ 1         │ 1           │ Claude et al. 1983 │
│ 3   │ γ      │ 6.0     │ 3.0           │ mol      │ Normal{Float64}(μ=1.0, σ=2.0)    │ The # of moles │ (-1, 1)     │ 0         │ 1           │ Dusse et al. 2000  │
```

Note that there is no check that the metadata you give is consistent.
These metadata will hopefully be useful for advanced usage of AIBECS, e.g., using prior information and/or bounds for optimization.
"""
abstract type AbstractParameters{T} <: AbstractVector{T} end
# TODO: Should it still be a subtype of AVec?
const APar = AbstractParameters # alias to make code more readable

"""
    symbols(p)

Returns the symbols in `p`.

Can also be used directly on the type of `p`
(because the symbols of `p::T` are contained in the type `T`).
"""
symbols(::Type{T}) where {T <: APar} = fieldnames(T)
symbols(::T) where {T <: APar} = symbols(T)

"""
    flattenable_symbols(p)

Returns the flattenable symbols in `p`.

The flattenable symbols are those symbols that are kepth when using `p` as a vector,
e.g., when doing `vec(p)`.
(Useful when passing parameters to an optimization routine expeting a vector of optimizable parameters.)

Can also be used directly on the type of `p`
(because the flattenable symbols of `p::T` are contained in the type `T`).
"""
flattenable_symbols(::T) where {T <: APar} = flattenable_symbols(T)
flattenable_symbols(::Type{T}) where {T <: APar} = fieldnames(T)[collect(flattenable(T))]

"""
    values(p::T) where {T <: AbstractParameters}

Returns a vector of **all** the values of `p`.

Note that `values(p)` is different from `vec(p)`.
"""
Base.values(p::T) where {T <: APar} = [getfield(p, s) for s in symbols(T)]

"""
    flattenable_values(p::T) where {T <: AbstractParameters}

Returns a vector of the **flattenable** values of `p`.

Note that `vec(p)` is different from `values(p)`.
"""
flattenable_values(p::T) where {T <: APar} = [getfield(p, s) for s in flattenable_symbols(T)]

"""
    vec(p::T) where {T <: AbstractParameters}

Returns a **SI-unit-converted** vector of flattenable values of `p`.

Note that `vec(p) ≠ flattenable_values(p)` if `p` has units.
"""
Base.vec(p::T) where {T <: APar} = [UnPack.unpack(p, Val(s)) for s in flattenable_symbols(T)]


"""
    length(p::AbstractParameter)

Returns the length of the **flattened/optimzable** vector of `p`.

May be different from the number of parameters.
Can also be used directly on the type of `p`.
"""
Base.length(::T) where {T <: APar} = sum(flattenable(T))
Base.length(::Type{T}) where {T <: APar} = sum(flattenable(T))

"""
    size(p::AbstractParameter)

Returns the size of the **flattened/optimzable** vector of `p`.

May be different from the number of parameters.
Can also be used directly on the type of `p`.
"""
Base.size(::T) where {T <: APar} = (length(T),)
Base.size(::Type{T}) where {T <: APar} = (length(T),)

function Base.show(io::IO, p::T) where T <: APar
    print(T)
    t = table(p)
    show(io, t; summary=false)
end
function Base.show(io::IO, m::MIME"text/plain", p::T) where T <: APar
    print(T)
    t = table(p)
    show(io, m, t; summary=false)
end

function table(p::T, nf::NamedTuple) where {T <: APar}
    t = DataFrame(Symbol = collect(symbols(p)), Value = values(p))
    for (n,f) in zip(keys(nf), nf)
        setproperty!(t, n, collect(f(p)))
    end
    return t
end

"""
    table(p)

Returns a `DataFrame` (a table) of `p`.
Useful for printing and saving into an actual text/latex table.
"""
function table(p::APar)
    ks, vs = Symbol[], Function[]
    all(isnothing, initial_value(p)) || (push!(ks, Symbol("Initial value")); push!(vs, initial_value))
    all(isequal(1), units(p)) || (push!(ks, :Unit); push!(vs, units))
    all(isnothing, prior(p)) || (push!(ks, Symbol("Prior")); push!(vs, prior))
    all(isnothing, limits(p)) || (push!(ks, Symbol("Limits")); push!(vs, limits))
    all(isempty, description(p)) || (push!(ks, Symbol("Description")); push!(vs, description))
    all(isnothing, bounds(p)) || (push!(ks, Symbol("Bounds")); push!(vs, bounds))
    any(logscaled(p)) && (push!(ks, Symbol("Logscaled")); push!(vs, logscaled))
    all(flattenable(p)) || (push!(ks, Symbol("Optimizable")); push!(vs, flattenable))
    all(isnothing, reference(p)) || (push!(ks, Symbol("Reference")); push!(vs, reference))
    nf = (; zip(ks, vs)...)
    t = table(p, nf)
end


"""
    latex(p)

Returns a LaTeX-formatted table of the parameters.
"""
latex(p::APar) = show(stdout, MIME("text/latex"), table(p))

"""
    unpack(p <: AbstractParameters, s)

Unpacks the parameter `s` from `p`.

Note this is specialized and will convert the parameter value to SI units.
"""
@inline UnPack.unpack(p::T, ::Val{f}) where {T<:APar,f} = ustrip(upreferred(getproperty(p, f) * units(T, f)))

"""
    +(p::T, v::Vector) where {T <: AbstractParameters}

Adds the flattened vector `v` to `p`.

**Warning:** This method for `+` is implemented only for differentiation
using dual and hyperdual numbers.
If you want to change the values of `p`, you should do so explicitly
rather than use this `+` method.
"""
Base.:+(p::T, v::Vector) where {T <: APar} = reconstruct(T, flattenable_values(p) + v)

"""
    reconstruct(T, v)

Reconstructs the parameter of type `T` from flattenable vector `v`.
"""
function reconstruct(::Type{T}, v::Vector{V}) where {T <: APar, V}
    all(isnothing, initial_value(T)) && error("Can't reconstruct without initial values")
    reconstructed_v = convert(Vector{V}, collect(initial_value(T))) # fill in initial values
    reconstructed_v[collect(flattenable(T))] .= v            # fill in v
    return T.name.wrapper{V}(reconstructed_v...)  # wrapper allows to change type
end

"""
    getindex(p::T, i) where {T <: AbstractParameters}

Returns the i-th element of vec(p).

This is not efficient and only used for testing the derivatives with ForwardDiff.
"""
Base.getindex(p::T, i) where {T <: APar} = flattenable_values(p)[i]


function (::Type{T})(args::Quantity...) where {T <: APar}
    all(isequal(1), units(T)) && error("$T needs `units` for this construction to work")
    return T([ustrip(x |> units(T, f)) for (x,f) in zip(args, fieldnames(T))]...)
end

function (::Type{T})(;kwargs...) where {T <: APar}
    all(isnothing, initial_value(T)) && length(kwargs) ≠ length(symbols(T)) && error("$T needs `initial_value` if one of the parameters is not supplied")
    value(f::Symbol, v::Quantity) = ustrip(v |> units(T, f))
    value(f::Symbol, v) = v
    return T([f ∈ keys(kwargs) ? value(f, kwargs[f]) : initial_value(T, f) for f in fieldnames(T)]...)
end

#===============================
Writing to savable formats
===============================#

function Base.Dict(p::T, s=symbols(p)) where {T<:APar}
    v = [getfield(p,s) for s in s]
    u = [units(p,s) for s in s]
    return Dict([(s,v*u) for (s,v,u) in zip(s,v,u)])
end
function Base.NamedTuple(p::T, s=symbols(p)) where {T<:APar}
    v = [getfield(p,s) for s in s]
    u = [units(p,s) for s in s]
    return (; zip(s, v .* u)...)
end

#=====================
mismatch of parameters
=====================#

# in p space (should not be needed except for checking/testing)
mismatch(p::T, k) where {T<:APar} = mismatch(T, bijector(T, k)(p))
mismatch(p::T) where {T<:APar} = mismatch(T, bijector(T)(p))

# in λ space
mismatch(::Type{T}, λ) where {T<:APar} = sum(mismacth(T, kᵢ, λᵢ) for (k, λᵢ) in zip(flattenable_symbols(T), λ))
mismatch(::Type{T}, k, λ) where {T<:APar} = mismatch(prior(T, k), λ)
mismatch(d::ContinuousUnivariateDistribution, λ) = -logpdf(transformed(d), λ)

# TODO check if I need that below:
# The functions below is just for ForwardDiff to work with vectors instead of p
# which requires `mismatch` to know about the priors, which are containted in `T`
function generate_objective(ωs, μx, σ²x, v, ωp, ::Type{T}) where {T<:APar}
    nt, nb = length(ωs), length(v)
    tracers(x) = state_to_tracers(x, nb, nt)
    f(x, p) = ωp * mismatch(T, p) +
        sum([ωⱼ * mismatch(xⱼ, μⱼ, σⱼ², v) for (ωⱼ, xⱼ, μⱼ, σⱼ²) in zip(ωs, tracers(x), μx, σ²x)])
    return f
end

#==================
Change of variables
==================#


# using Bijectors and evaluating parameter mismatch in λ-space
# Usage should look like this
# F, ∇ₓF = generate_state_function(...)
# f, ∇ₓf = generate_objective(...) # This requires to dispatch them like f(p::APar) = ... and f(λ::Vector) = ...
# obj, grad, hess = F1Method(F, ∇ₓF, f, ∇ₓf, ...) # directly in λ-space
# p2λ = bijector(MyParams) # just reexport bijector and inverse
# λ2p = inverse(p2λ)
# λ₀ = p2λ(p₀)
# results = optimize(obj, grad, hess, λ₀, ...)
# λ_opt = results.minimizer
# p_opt = λ2p(λ_opt) # and done and dusted...
# some actual code that I should write and test
import Bijectors: bijector
"""
    bijector(T::AbstractParameters [, k::Symbol])

Returns the function for the change of variables of the parameters.

The function is a bijection from the supports/domains of the priors to ℝⁿ,
from the Bijectors.jl package.
You can specify the parameter symbol to get the bijector of that parameter.
"""
bijector(::Type{T}) where {T<:APar} = bijector(prior(T)[flattenable(T)]) # should work out of the box!
bijector(::Type{T}, k) where {T<:APar} = bijector(prior(T, k)) # not sure even needed
export bijector

"""
p2λ(p::AbstractParameters)

Converts `p` to a real-valued vector for optimization.
(referred to as the λ space in AIBECS)
"""
p2λ(p::T) = bijector(T)(flattenable_values(p))
"""
λ2p(T::Type{AbstractParameters}, λ::Vector)

Converts real-valued vector `λ` back to parameters object `p`.

Note that the instance of your parameters type `T` is required here because
it contains information on non-optimizable parameters and priors of optimizable parameters
"""
λ2p(T::Type{APar}, λ) = T(; zip(flattenable_symbols(T), inverse(bijector(T))(λ))...)
export p2λ, λ2p

export AbstractParameters, latex

