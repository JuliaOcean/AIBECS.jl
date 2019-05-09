# Set of convenient tools to convert a table of parameters into the right type and
# the default parameter

#=============================================
Example use:
t = empty_parameter_table()    # initialize table of parameters
add_parameter!(t, :umax, 1e-3u"mmol/m^3/d",
    optimizable = true,
    description = "Maximum uptake rate (Michaelis-Menten)",
    LaTeX = "U_\\mathrm{max}")
add_parameter!(t, :ku, 1e-1u"mmol/m^3",
    optimizable = true,
    description = "Half-saturation constant (Michaelis-Menten)",
    LaTeX = "k_\\vec{u}")
add_parameter!(t, :w₀, 1u"m/d",
    optimizable = true,
    description = "Sinking velocity at surface",
    LaTeX = "w_0")
add_parameter!(t, :w′, 1u"d^-1",
    optimizable = true,
    description = "Vertical gradient of sinking velocity",
    LaTeX = "w'")
add_parameter!(t, :κ, 0.25u"d^-1",
    optimizable = true,
    description = "Remineralization rate constant",
    LaTeX = "\\kappa")
add_parameter!(t, :DINgeo, 2.3u"mmol/m^3",
    optimizable = true,
    description = "Mean DIN concentration",
    LaTeX = "\\xNUT^\\obs")
add_parameter!(t, :τg, 1u"Myr",
    description = "Geological restoring timescale",
    LaTeX = "\\tau_\\mathrm{geo}")
add_parameter!(t, :ω, 1e-4u"1",
    description = "Relative weight of cost of parameters (Hyper parameter)",
    LaTeX = "\\omega")
initialize_parameter_type(t)   # Generate the parameter type
=============================================#

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
Generate Para
====#

import Flatten: flattenable

macro make_struct(struct_name, schema...)
    fields = [:($(x[1])::U | $(x[2])) for x in schema...]
    esc(quote @flattenable mutable struct $struct_name{U} <: AbstractVector{U}
        $(fields...)
        end
    end)
end


function initialize_parameter_type(t)
    symbols = t[:symbol]
    optimizables = t[:optimizable]
    optsymbols = symbols[optimizables]
    schema = [(Symbol(x),Symbol(y)) for (x,y) in zip(symbols, optimizables)]
    eval( :(@make_struct $(:Para) $schema))

    printunits = t[:printunit]
    baseunits = t[:unit]
    values = t[:value]

    μs = [μ for (μ, opt) in zip(t[:mean_obs], optimizables) if opt]
    σ²s = [σ² for (σ², opt) in zip(t[:variance_obs], optimizables) if opt]

    @eval begin
        # Printing functionality
        function Base.show(io::IO, p::Para)
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
        Base.show(io::IO, ::MIME"text/plain", p::Para) = Base.show(io, p)
        # constants
        const p₀ = Para($values...)
        const m = length(fieldnameflatten(p₀))
        const m_all = length(fieldnames(typeof(p₀)))
        const mean_pobs = $μs
        const variance_pobs = $σ²s
        export p₀, m, mean_pobs, variance_pobs, m_all
        # overloads
        Base.length(p::Para) = length(fieldnameflatten(p))
        # Make Para an iterable for the parameters to be able to `collect` it into a vector
        Base.iterate(p::Para, i=1) = i > m_all ? nothing : (getfield(p, i), i + 1)
        # Convert p to a vector and vice versa
        Base.vec(p::Para) = collect((p...,))
        Base.copy(p::Para) = Para(vec(p)...)
        Base.convert(T, p) = Para(convert(Vector{T}, vec(p))...)
        opt_para(p, v) = Flatten.reconstruct(p, v)
        function opt_para(p::Para{Tₚ}, v::Vector{Tᵥ}) where {Tₚ, Tᵥ}
            Flatten.reconstruct(convert(Tᵥ, p), v)
        end
        opt_para(v) = opt_para(p₀, v)
        optvec(p::Para) = flatten(Vector, p)
        export optvec
        # Testing equality and approx
        Base.:≈(p₁::Para, p₂::Para) = vec(p₁) ≈ vec(p₂)
        Base.:(==)(p₁::Para, p₂::Para) = vec(p₁) == vec(p₂)
        # Overloads for being a subtype of Vector
        strerror = "Index of of bounds!"
        Base.getindex(p::Para, i::Int) = i < 1 || i > m ? error(strerror) : getfield(p, $optsymbols[i])
        Base.setindex!(p::Para, v, i::Int) = i < 1 || i > m ? error(strerror) : setfield!(p, $optsymbols[i], v)
        # base overloads
        Base.:+(p::Para, v::Vector) = opt_para(p₀, optvec(p) + v)
        Base.:-(p::Para, v::Vector) = opt_para(p₀, optvec(p) - v)
        Base.:+(p₁::Para, p₂::Para) = opt_para(p₀, optvec(p₁) + optvec(p₂))
        Base.:-(p₁::Para, p₂::Para) = opt_para(p₀, optvec(p₁) - optvec(p₂))
        Base.:*(s::Number, p::Para) = opt_para(p₀, s * optvec(p))
        Base.:*(p::Para, s::Number) = s * p
    end
end
export initialize_parameter_type

#function make_Para(t, v)
#    make_Para_struct(t)
#    return Para(v...)
#end

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


