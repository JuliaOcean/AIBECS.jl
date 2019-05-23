
function empty_parameter_table()
    return DataFrame(
        symbol = Symbol[],
        value = Float64[],
        unit = Unitful.Units[],
        printunit = Unitful.Units[],
        mean_obs = Float64[],
        variance_obs = Float64[],
        optimizable = Bool[],
        description = String[],
        LaTeX = String[]
    )
end
export empty_parameter_table

"""
    new_parameter(symbol::Symbol,
                  quantity;
                  mean_obs=ustrip(upreferred(quantity)),
                  variance_obs=ustrip(upreferred(quantity))^2,
                  optimizable=false,
                  description="",
                  LaTeX="")

Creats a parameter (to be added to the parameters table).
If keyword argument `optimizable = false`, then observation mean and
variance are set to `NaN`.
Otherwise, these are set to keyword arguments `mean_obs` (and `variance_obs`)
if supplied, or to `quantity` (and its square), after converting it to
the preferred unit and stripping it of said unit if not.
Example: TODO
"""
new_parameter(symbol::Symbol,
              quantity;
              mean_obs=ustrip(upreferred(quantity)),
              variance_obs=ustrip(upreferred(quantity))^2,
              optimizable=false,
              description="",
              LaTeX="") = [symbol,
                           ustrip(upreferred(quantity)),
                           unit(upreferred(quantity)),
                           unit(quantity),
                           optimizable ? mean_obs : NaN,
                           optimizable ? variance_obs : NaN,
                           optimizable,
                           description,
                           LaTeX]
export new_parameter

"""
    add_parameter!(t::DataFrame, args...; kwargs...)

Adds a parameter to the parameters table `t`.
If keyword argument `optimizable = false`, then observation mean and
variance are set to `NaN`.
Otherwise, these are set to keyword arguments `mean_obs` (and `variance_obs`)
if supplied, or to `quantity` (and its square), after converting it to
the preferred unit and stripping it of said unit if not.
Example: TODO
Note for future edit of the docs: Don't repeat yourself between add and new param functions
"""
function add_parameter!(t::DataFrame, args...; kwargs...)
    if any(t[:symbol] .== args[1])
        error("Parameter $(args[1]) already exists! (Maybe delete it first?)")
    else
        push!(t, new_parameter(args...; kwargs...))
    end
end
export add_parameter!

delete_parameter!(t::DataFrame, i) = deleterows!(t, i)
function delete_parameter!(t::DataFrame, s::Symbol)
    i = findfirst(t[:symbol] .== s)
    if isnothing(i)
        println("Parameter $s did not exist. No deletions")
        return t
    else
        deleterows!(t, i)
    end
end
export delete_parameter!

#====
Generate Parameters
====#

import Flatten: flattenable

macro make_struct(struct_name, schema...)
    fields = [:($(x[1])::U | $(x[2])) for x in schema...]
    esc(quote @flattenable mutable struct $struct_name{U} <: AbstractVector{U}
        $(fields...)
        end
    end)
end


function initialize_Parameters_type(t)
    symbols = t[:symbol]
    optimizables = t[:optimizable]
    optsymbols = symbols[optimizables]
    schema = [(Symbol(x),Symbol(y)) for (x,y) in zip(symbols, optimizables)]
    eval( :(@make_struct $(:Parameters) $schema))

    printunits = t[:printunit]
    baseunits = t[:unit]
    values = t[:value]

    μs = [μ for (μ, opt) in zip(t[:mean_obs], optimizables) if opt]
    σ²s = [σ² for (σ², opt) in zip(t[:variance_obs], optimizables) if opt]

    @eval begin
        # Printing functionality
        function Base.show(io::IO, p::Parameters)
            println(typeof(p))
            for (s, pu, bu, opt) in zip($symbols, $printunits, $baseunits, $optimizables)
                v = getfield(p, s)
                if pu == 1
                    val, ppu = v, ""
                else
                    val, ppu = ustrip(uconvert(pu, v * bu)), unicodify(pu)
                end
                optstr = opt ? "" : "(fixed)"
                print_type(io, s, val, ppu, optstr)
            end
        end
        Base.show(io::IO, ::MIME"text/plain", p::Parameters) = Base.show(io, p)
        # constants
        const p₀ = Parameters($values...)
        const m = length(fieldnameflatten(p₀))
        const m_all = length(fieldnames(typeof(p₀)))
        const mean_pobs = $μs
        const variance_pobs = $σ²s
        export p₀, m, mean_pobs, variance_pobs, m_all
        # overloads
        Base.length(p::Parameters) = length(fieldnameflatten(p))
        Base.size(p::Parameters) = (length(p),) # ForwardDiff requirement
        # Make Parameters an iterable for the parameters to be able to `collect` it into a vector
        Base.iterate(p::Parameters, i=1) = i > m_all ? nothing : (getfield(p, i), i + 1)
        # Convert p to a vector and vice versa
        Base.vec(p::Parameters) = collect((p...,))
        Base.copy(p::Parameters) = Parameters(vec(p)...)
        Base.convert(::Type{Parameters{T1}}, p::Parameters{T2}) where {T1, T2} = Parameters(convert(Vector{T1}, vec(p))...)
        Base.convert(::Type{Parameters{T}}, p::Parameters{T}) where T = p
        opt_para(p, v) = Flatten.reconstruct(p, v)
        function opt_para(p::Parameters{Tₚ}, v::Vector{Tᵥ}) where {Tₚ, Tᵥ}
            Flatten.reconstruct(convert(Parameters{Tᵥ}, p), v)
        end
        opt_para(v) = opt_para(p₀, v)
        optvec(p::Parameters) = flatten(Vector, p)
        optvec(v) = v # ForwardDiff requirement
        export optvec
        # Testing equality and approx
        Base.:≈(p₁::Parameters, p₂::Parameters) = vec(p₁) ≈ vec(p₂)
        Base.:(==)(p₁::Parameters, p₂::Parameters) = vec(p₁) == vec(p₂)
        # Overloads for being a subtype of Vector
        strerror = "Index of of bounds!"
        Base.getindex(p::Parameters, i::Int) = i < 1 || i > m ? error(strerror) : getfield(p, $optsymbols[i])
        Base.setindex!(p::Parameters, v, i::Int) = i < 1 || i > m ? error(strerror) : setfield!(p, $optsymbols[i], v)
        # base overloads
        Base.:+(p::Parameters, v::Vector) = opt_para(p₀, optvec(p) + v)
        Base.:-(p::Parameters, v::Vector) = opt_para(p₀, optvec(p) - v)
        Base.:+(p₁::Parameters, p₂::Parameters) = opt_para(p₀, optvec(p₁) + optvec(p₂))
        Base.:-(p₁::Parameters, p₂::Parameters) = opt_para(p₀, optvec(p₁) - optvec(p₂))
        Base.:*(s::Number, p::Parameters) = opt_para(p₀, s * optvec(p))
        Base.:*(p::Parameters, s::Number) = s * p
    end
end
export initialize_Parameters_type


#=======
Printing
=======#

print_type(io, f, val::Float64, ppu, s) = @printf io "%6s = %8.2e [%s] %s\n" f val ppu s
print_type(io, f, val::Dual{Float64}, ppu, s) = @printf io "%6s = %8.2e + %8.2eε [%s] %s\n" f ℜ(val) 𝔇(val) ppu s
print_type(io, f, val::Complex{Float64}, ppu, s) = @printf io "%6s = %8.2e + %8.2ei [%s] %s\n" f ℜ(val) ℑ(val) ppu s
print_type(io, f, val::Hyper{Float64}, ppu, s) = @printf io "%6s = %8.2e + %8.2eε₁ + %8.2eε₂ + %8.2eε₁ε₂ [%s] %s\n" f ℜ(val) ℌ₁(val) ℌ₂(val) ℌ(val) ppu s
ℜ(x::Complex) = real(x)
ℑ(x::Complex) = imag(x)
ℜ(x::Dual) = DualNumbers.realpart(x)
𝔇(x::Dual) = DualNumbers.dualpart(x)
ℜ(x::Hyper) = HyperDualNumbers.realpart(x)
ℌ(x::Hyper) = HyperDualNumbers.ε₁ε₂part(x)
ℌ₁(x::Hyper) = HyperDualNumbers.ε₁part(x)
ℌ₂(x::Hyper) = HyperDualNumbers.ε₂part(x)

function unicodify(U::Unitful.Units)
    str = string(U)
    str = replace(str, r"\^-1" => s"⁻¹")
    str = replace(str, r"\^-2" => s"⁻²")
    str = replace(str, r"\^-3" => s"⁻³")
    str = replace(str, r"\^1" => s"¹")
    str = replace(str, r"\^2" => s"²")
    str = replace(str, r"\^3" => s"³")
    return str
end


