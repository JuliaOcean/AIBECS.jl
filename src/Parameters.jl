# Set of convenient tools to convert a table of parameters into the right type and
# the default parameter

#=============================================
Format of DataFrame to be supplied by the user

julia> Parameters_table = DataFrame()

julia> Parameters_table.symbol =
    [
)
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

function new_parameter(
        symbol::Symbol,
        quantity::Quantity;
        mean_obs=NaN,
        variance_obs=NaN,
        optimizable=false,
        description="",
        LaTeX=""
    )
    if optimizable && (isnan(mean_obs) || isnan(variance_obs))
        error("If the parameter is optimizable you need to give a prior mean and variance!")
    else
        return [symbol,
            ustrip(upreferred(quantity)),
            unit(upreferred(quantity)),
            unit(quantity),
            mean_obs,
            variance_obs,
            optimizable,
            description,
            LaTeX]
    end
end
export new_parameter

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
Dummy generation
=====#
#t = empty_parameter_table()
#add_parameter!(t, :α, 1.0u"mmol/m^3")
#add_parameter!(t, :β, 2.0u"mmol/m^3"; mean_obs = 2.0, variance_obs = 3.0, optimizable = false)
#add_parameter!(t, :ϵ, 3.0u"mmol/m^3"; mean_obs = 2.0, variance_obs = 3.0, optimizable = true)
#add_parameter!(t, :τ, 4.0u"mmol/l"; mean_obs = 2.0, variance_obs = 3.0, optimizable = false)
#add_parameter!(t, :γ, 5.0u"mol/m^3"; mean_obs = 2.0, variance_obs = 3.0, optimizable = true)
#add_parameter!(t, :λ, 6.0u"mmol/m^3/s"; mean_obs = 2.0, variance_obs = 3.0, optimizable = false)
#add_parameter!(t, :Γ, 7.0u"mmol/m^3/d"; mean_obs = 2.0, variance_obs = 3.0, optimizable = true)
#add_parameter!(t, :Ω, 8.0u"mmol/m^2"; mean_obs = 2.0, variance_obs = 3.0, optimizable = true)
#====
end Dummy generation
=====#


#====
Generate Constants
====#

#function generate_constants(t)
#    for i in findall(t[:optimizable] .== false)
#        @eval $(t[:symbol][i]) = $(t[:value][i])
#    end
#end
#
#macro generate_Para(t)
#    filter(row -> row[:optimizable] == true, t)[:value]
#end
#
#macro make_Para(t)
#    println("Getting there?")
#    topt = $t[:optimizable]
#    println("Maybe?")
#    iopt = findall(topt)
#    println(iopt)
#    fields = [:( $(t[:symbol][i])::U ) for i in iopt]
#    println(fields)
#    esc(quote @with_kw struct Para{U} <: AbstractVector{U}
#        $(fields...)
#        end
#    end)
#end

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
        export p₀
        const m = length(fieldnameflatten(p₀))
        export m
        const m_all = length(fieldnames(typeof(p₀)))
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
        # Testing equality and approx
        Base.:≈(p₁::Para, p₂::Para) = vec(p₁) ≈ vec(p₂)
        Base.:(==)(p₁::Para, p₂::Para) = vec(p₁) == vec(p₂)
        # Overloads for being a subtype of Vector
        strerror = "Index of of bounds!"
        Base.getindex(p::Para, i::Int) = i < 1 || i > m ? error(strerror) : getfield(p, $optsymbols[i])
        Base.setindex!(p::Para, v, i::Int) = i < 1 || i > m ? error(strerror) : setfield!(p, $optsymbols[i], v)
        # Convert p to λ and vice versa, needed by AIBECS!
        optvec(p::Para) = flatten(Vector, p)
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

#====
Dummy generation
=====#
#make_Para_struct(t)
#p₀
#v1 = vec(p₀)
#v2 = optvec(p₀)
#p₀ + v2
#p₀ - v2 / 2
#p₀ + ε * ones(4)
#3p₀
#ε * p₀
#p₀[1]
#p₀[4]
#p₀ + ε₁ * v2
#====
end Dummy generation
=====#
