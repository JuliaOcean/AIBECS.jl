using Parameters, FieldMetadata, Flatten, DataFrames, Unitful, Distributions
import Flatten: flattenable
@metadata initial_value nothing
import FieldMetadata: @units, @reunits, units
import FieldMetadata: @prior, @reprior, prior
import FieldMetadata: @description, @redescription, description
@metadata bounds nothing
import FieldMetadata: @logscaled, @relogscaled, logscaled
import FieldMetadata: @flattenable, @reflattenable, flattenable
@eval const $(Symbol("@optimizable")) = $(Symbol("@flattenable"))
optimizable = flattenable
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

- nice syntax for unpacking parameters in functions via Parameters' `@unpack` macro
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
julia> @with_kw struct SimpleParams{T} <: AbstractParameters{T}
           α::T
           β::T
           γ::T
       end
SimpleParams
```

To create an instance of the `SimpleParams(Float64)` type, you can do

```jldoctest
julia> SimpleParams(α=1.0, β=2.0, γ=3.0)
SimpleParams{Float64}
│ Row │ Symbol │ Value   │
│     │ Symbol │ Float64 │
├─────┼────────┼─────────┤
│ 1   │ α      │ 1.0     │
│ 2   │ β      │ 2.0     │
│ 3   │ γ      │ 3.0     │
```

The parameter types you define in AIBECS can contain default values, if you want then to, 
thanks to the Parameters.jl interface. For example, 

```jldoctest
julia> @with_kw struct SimpleParamsWithDefault{T} <: AbstractParameters{T}
           α::T = 1.0
           β::T = 2.0
           γ::T = 3.0
       end
SimpleParamsWithDefault
```

This allows for very simple generation of instances of that type.

```jldoctest
julia> SimpleParamsWithDefault()
SimpleParamsWithDefault{Float64}
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

julia> simplef(SimpleParamsWithDefault()) # 1.0 + 3.0
4.0
```

More complex examples are permitted by adding metadata (thanks to FieldMetadata.jl).
You can add units

```jldoctest
julia> @units @with_kw struct UnitParams{T} <: AbstractParameters{T}
           α::T = 1.0 | u"km"
           β::T = 2.0 | u"hr"
           γ::T = 3.0 | u"mol"
       end ;

julia> UnitParams()
UnitParams{Float64}
│ Row │ Symbol │ Value   │ Unit     │
│     │ Symbol │ Float64 │ Unitful… │
├─────┼────────┼─────────┼──────────┤
│ 1   │ α      │ 1.0     │ km       │
│ 2   │ β      │ 2.0     │ hr       │
│ 3   │ γ      │ 3.0     │ mol      │
```

Note that when adding units to your parameters, they will be converted to SI
when unpacked, as in, e.g.,

```jldoctest
julia> function speed(p::UnitParams)
           @unpack α, β = p
           return α / β
       end
speed (generic function with 1 method)

julia> speed(UnitParams()) # 1.0 km / 2.0 hr in m/s
0.1388888888888889
```

Thanks to the FieldMetaData interface, you can chain the following preloaded metadata:

- initial_value
- units (from Unitful.jl and UnitfulAstro.jl)
- prior (from Distributions.jl)
- description (`String`)
- bounds (2-element `Tuple`)
- logscaled (`Bool`)
- optimizable (alias for `flattenable`)
- reference (`String`)

Here is an example of parameter with all the possible metadata available in AIBECS:

```jldoctest
julia> @reference @optimizable @logscaled @bounds @description @prior @units @initial_value @with_kw struct FullParams{T} <: AbstractParameters{T}
           α::T = 1.0 | 0.0 | u"km"  | Normal(0,1)    | "The distance"   | (-Inf, Inf) | false | false | "Jean et al., 2042" 
           β::T = 2.0 | 0.0 | u"hr"  | LogNormal(0,1) | "The time"       | (   0, Inf) | true  | true  | "Claude et al. 1983" 
           γ::T = 3.0 | 0.0 | u"mol" | Normal(1,2)    | "The # of moles" | (  -1,   1) | false | true  | "Dusse et al. 2000"
       end ;

julia> FullParams()
FullParams{Float64}
│ Row │ Symbol │ Value   │ Initial value │ Unit     │ Prior                            │ Description    │ Bounds      │ Logscaled │ Optimizable │ Reference          │
│     │ Symbol │ Float64 │ Float64       │ Unitful… │ Distribu…                        │ String         │ Tuple…      │ Bool      │ Bool        │ String             │
├─────┼────────┼─────────┼───────────────┼──────────┼──────────────────────────────────┼────────────────┼─────────────┼───────────┼─────────────┼────────────────────┤
│ 1   │ α      │ 1.0     │ 0.0           │ km       │ Normal{Float64}(μ=0.0, σ=1.0)    │ The distance   │ (-Inf, Inf) │ 0         │ 0           │ Jean et al., 2042  │
│ 2   │ β      │ 2.0     │ 0.0           │ hr       │ LogNormal{Float64}(μ=0.0, σ=1.0) │ The time       │ (0, Inf)    │ 1         │ 1           │ Claude et al. 1983 │
│ 3   │ γ      │ 3.0     │ 0.0           │ mol      │ Normal{Float64}(μ=1.0, σ=2.0)    │ The # of moles │ (-1, 1)     │ 0         │ 1           │ Dusse et al. 2000  │
```

Note that there is no check that the metadata you give is consistent.
These metadata will hopefully be useful for advanced usage of AIBECS, e.g., using prior information and/or bounds for optimization.
"""
abstract type AbstractParameters{T} <: AbstractVector{T} end

"""
    symbols(p)

Returns the symbols in `p`.

Can also be used directly on the type of `p`
(because the symbols of `p::T` are contained in the type `T`).
"""
symbols(::Type{T}) where {T <: AbstractParameters} = fieldnames(T)
symbols(::T) where {T <: AbstractParameters} = fieldnames(T)

"""
    flattenable_symbols(p)

Returns the flattenable symbols in `p`.

The flattenable symbols are those symbols that are kepth when using `p` as a vector,
e.g., when doing `vec(p)`.
(Useful when passing parameters to an optimization routine expeting a vector of optimizable parameters.)

Can also be used directly on the type of `p`
(because the flattenable symbols of `p::T` are contained in the type `T`).
"""
flattenable_symbols(::T) where {T <: AbstractParameters} = fieldnames(T)[collect(flattenable(T))]
flattenable_symbols(::Type{T}) where {T <: AbstractParameters} = fieldnames(T)[collect(flattenable(T))]

"""
    values(p)

Returns a vector of **all** the values of `p`.

Note that `values(p)` is different from `vec(p)`.
"""
values(p::T) where {T <: AbstractParameters} = [getfield(p, s) for s in symbols(T)]

"""
    flattenable_values(p)

Returns a vector of the **flattenable** values of `p`.

Note that `vec(p)` is different from `values(p)`.
"""
flattenable_values(p::T) where {T <: AbstractParameters} = [getfield(p, s) for s in flattenable_symbols(T)]

"""
    vec(p)

same as `flattenable_values(p)`.
"""
vec(p::AbstractParameters) = flattenable_values(p)

"""
    length(p::AbstractParameter)

Returns the length of the **flattened** vector of `p`.

Can also be used directly on the type of `p`.
"""
Base.length(::T) where {T <: AbstractParameters} = sum(flattenable(T))
Base.length(::Type{T}) where {T <: AbstractParameters} = sum(flattenable(T))

"""
    size(p::AbstractParameter)

Returns the size of the **flattened** vector of `p`.

Can also be used directly on the type of `p`.
"""
Base.size(::T) where {T <: AbstractParameters} = (length(T),)
Base.size(::Type{T}) where {T <: AbstractParameters} = (length(T),)

strerror = "Index out of bounds!"
Base.getindex(p::T, i::Int) where {T <: AbstractParameters} = Parameters.unpack(p, Val(flattenable_symbols(T)[i]))

function Base.show(io::IO, p::T) where T <: AbstractParameters
    print(T)
    t = DataFrame(symbol = collect(parameter_symbols(p)), value = vec(p))
    show(io, t; summary=false)
end
function Base.show(io::IO, m::MIME"text/plain", p::T) where T <: AbstractParameters
    print(T)
    t = table(p)
    show(io, m, t; summary=false)
end

function table(p::AbstractParameters, nf::NamedTuple)
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
function table(p::AbstractParameters)
    ks, vs = Symbol[], Function[]
    all(isnothing, initial_value(p)) || (push!(ks, Symbol("Initial value")); push!(vs, initial_value))
    all(isequal(1), units(p)) || (push!(ks, :Unit); push!(vs, units))
    all(isnothing, prior(p)) || (push!(ks, Symbol("Prior")); push!(vs, prior))
    all(isempty, description(p)) || (push!(ks, Symbol("Description")); push!(vs, description))
    all(isnothing, bounds(p)) || (push!(ks, Symbol("Bounds")); push!(vs, bounds))
    any(logscaled(p)) && (push!(ks, Symbol("Logscaled")); push!(vs, logscaled))
    all(optimizable(p)) || (push!(ks, Symbol("Optimizable")); push!(vs, optimizable))
    all(isnothing, reference(p)) || (push!(ks, Symbol("Reference")); push!(vs, reference))
    nf = (; zip(ks, vs)...)
    t = table(p, nf)
end


"""
    latex(p)

Returns a LaTeX-formatted table of the parameters.
"""
latex(p::AbstractParameters) = show(stdout, MIME("text/latex"), table(p))

#"""
#    getproperty(p::AbstractParameters, s)
#
#Returns the value of the parameter of symbol `s` in `p`.
#
#Note that the value returned is unitless, but if `p` has units, the value is converted to SI units beforehand.
#"""
#Base.getproperty(p::T, s::Symbol) where {T <: AbstractParameters} = ustrip(upreferred(getfield(p, s) * units(T, s)))

@inline Parameters.unpack(p::T, ::Val{f}) where {T<:AbstractParameters,f} = ustrip(upreferred(getproperty(p, f) * units(T, f)))


#=


=#


