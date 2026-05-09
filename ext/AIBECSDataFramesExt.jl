module AIBECSDataFramesExt

using AIBECS
using AIBECS: AbstractParameters, symbols, initial_value, units, prior,
    limits, description, bounds, logscaled, reference
using Flatten: flattenable
using DataFrames

function _table(p::T, nf::NamedTuple) where {T <: AbstractParameters}
    t = DataFrame(Symbol = collect(symbols(p)), Value = values(p))
    for (n, f) in zip(keys(nf), nf)
        setproperty!(t, n, collect(f(p)))
    end
    return t
end

function AIBECS.table(p::T) where {T <: AbstractParameters}
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
    return _table(p, nf)
end

AIBECS.latex(p::T) where {T <: AbstractParameters} = show(stdout, MIME("text/latex"), AIBECS.table(p))

end # module
