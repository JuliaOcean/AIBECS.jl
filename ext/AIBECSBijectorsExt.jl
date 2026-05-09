module AIBECSBijectorsExt

using AIBECS
using AIBECS: AbstractParameters, flattenable_symbols, flattenable_values, prior, reconstruct
import Bijectors: bijector
using Bijectors: inverse

bijector(::Type{T}) where {T <: AbstractParameters} = bijector.([p for p in prior(T) if !isnothing(p)])
bijector(::Type{T}, k) where {T <: AbstractParameters} = bijector(prior(T, k))

function AIBECS.p2λ(p::T) where {T <: AbstractParameters}
    return [bijector(T, k)(v) for (k, v) in zip(flattenable_symbols(p), flattenable_values(p))]
end

function AIBECS.λ2p(::Type{T}, λs) where {T <: AbstractParameters}
    ps = [inverse(bijector(T, k))(λ) for (k, λ) in zip(flattenable_symbols(T), λs)]
    return reconstruct(T, ps)
end

# parameter-space mismatch — needs `bijector` so it lives here
AIBECS.mismatch(p::T, k::Symbol) where {T <: AbstractParameters} = AIBECS.mismatch(T, bijector(T, k)(p))
AIBECS.mismatch(p::T) where {T <: AbstractParameters} = AIBECS.mismatch(T, AIBECS.p2λ(p))
AIBECS.mismatch(::Type{T}, p::Tp) where {T <: AbstractParameters, Tp <: AbstractParameters} = AIBECS.mismatch(p)

end # module
