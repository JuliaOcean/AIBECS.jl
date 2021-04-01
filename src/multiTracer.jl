



#= ============================================
Generate 𝐹 and ∇ₓ𝐹 from user input
============================================ =#

# TODO replace this with DiffEqOperators when possible
# SciMLBase
# ├─ AbstractSciMLOperator
# │  └─ AbstractDiffEqOperator
# │     ├─ AbstractDiffEqLinearOperator
# │     │  ├─ DiffEqIdentity
# │     │  ├─ DiffEqScalar
# │     │  ├─ DiffEqArrayOperator
# │     │  ├─ FactorizedDiffEqArrayOperator
# │     │  ├─ AbstractDerivativeOperator      <- DiffEqOperators
# │     │  │  └─ DerivativeOperator           <- DiffEqOperators
# │     │  ├─ AbstractDiffEqCompositeOperator <- DiffEqOperators
# │     │  │  ├─ DiffEqScaledOperator         <- DiffEqOperators
# │     │  │  ├─ DiffEqOperatorCombination    <- DiffEqOperators
# │     │  │  └─ DiffEqOperatorComposition    <- DiffEqOperators
# │     │  └─ AbstractMatrixFreeOperator      <- DiffEqOperators
# │     └─ AffineDiffEqOperator
# └─ AbstractDiffEqAffineOperator <- DiffEqOperators
# exported:
# - AffineDiffEqOperator
# - DiffEqScalar
# - DiffEqArrayOperator
# - DiffEqIdentity
#
# I think I need to simply use DiffEqArrayOperator and let composition magic happen
# thanks to AbstractDiffEqCompositeOperator.
import Base: size, convert
import SciMLBase: AbstractDiffEqOperator, update_coefficients!
struct BlockDiagDiffEqOperator{T,T1} <: AbstractDiffEqLinearOperator{T}
    ops::T1
    function BlockDiagDiffEqOperator{T}(ops) where T
        new{T,typeof(ops)}(ops)
    end
end
Base.size(L::BlockDiagDiffEqOperator) = size(L.ops[1]) .* length(L.ops)
function (L::BlockDiagDiffEqOperator)(u,p,t::Number)
    update_coefficients!(L,u,p,t)
    sum(A*u for A in L.ops)
end
function (L::BlockDiagDiffEqOperator)(du,u,p,t::Number)
    update_coefficients!(L,u,p,t)
    nt = length(L.ops)
    nb = size(L.ops[1], 1)
    for i in 1:nt
        idx = (i-1)*nb.+(1:nb)
        mul!(du[idx], L.ops[i], u[idx])
    end
end
function update_coefficients!(L::BlockDiagDiffEqOperator,u,p,t)
    for A in L.ops
        update_coefficients!(A,u,p,t)
    end
end
convert(::Type{AbstractMatrix}, L::BlockDiagDiffEqOperator) = blockdiag((convert(AbstractMatrix, A) for A in L.ops)...)


# If Parameters type if given, convert it into λ2p and pass that on
AIBECSFunction(Ts::Tuple, Gs::Tuple, ::Type{P}) where {P <: APar} = AIBECSFunction(Ts, Gs, λ -> λ2p(P, λ))
# Special case if only one tracer (no tuple) is given as arguments
AIBECSFunction(T::AbstractDiffEqLinearOperator, G::Function, λ2p=nothing) = AIBECSFunction((T,), (G,), λ2p)
# dispatch depending on number of args of G functions
function AIBECSFunction(Ts::Tuple, Gs::Tuple, λ2p=nothing)
    nt = length(Ts)
    nb = size(Ts[1], 1)
    tracers(u) = state_to_tracers(u, nb, nt)
    tracer(u, i) = state_to_tracer(u, nb, nt, i)
    if all(isinplace(G, nt + 2) for G in Gs) # +2 to account for dx and p
        AIBECSFunction_fromipGs(Ts, Gs, nt, nb, tracers, tracer, λ2p)
    else
        AIBECSFunction_fromoopGs(Ts, Gs, nt, nb, tracers, tracer, λ2p)
    end
end
# AIBECSFunction from out-of-place Gs
function AIBECSFunction_fromoopGs(Ts, Gs, nt, nb, tracers, tracer, λ2p)
    function G(du, u, p)
        du = copy(u)
        for j in 1:nt
            tracer(du, j) .= Gs[j](tracers(u)..., p)
        end
        du
    end
    return AIBECSFunction_from_ipG(Ts, G, nt, nb, tracers, tracer, λ2p)
end
# AIBECSFunction from in-place Gs
function AIBECSFunction_fromipGs(Ts, Gs, nt, nb, tracers, tracer, λ2p)
    function G(du, u, p)
        for j in 1:nt
            Gs[j](tracer(du, j), tracers(u)..., p)
        end
        du
    end
    return AIBECSFunction_from_ipG(Ts, G, nt, nb, tracers, tracer, λ2p)
end
# Given inplace big G(du, u, p) and these other functions, build the ODEFunction
# Note here G is the "big" G for all the tracers!
function AIBECSFunction_from_ipG(Ts, G, nt, nb, tracers, tracer, λ2p)
    function f(du, u, p, t)
        G(du, u, p)
        for j in 1:nt
            update_coefficients!(Ts[j], nothing, p, nothing) # necessary
            mul!(tracer(du, j), Ts[j], tracer(u, j), 1, -1) # du <- G(u) - T * u
        end
        du
    end
    f(du, u, λ::Vector, t) = f(du, u, λ2p(λ), t)
    f(u, p, t) = (du = copy(u); f(du, u, p, t); du)
    f(u, λ::Vector, t) = f(u, λ2p(λ), t)
    # i) Create ∇G(x) as a DiffEqArrayOperator
    jac0 = vcat((hcat((Diagonal(ones(nb)) for _ in 1:nt)...) for _ in 1:nt)...) # nt×nt blocks of diagonals
    colorsG = matrix_colors(jac0)
    # TODO, make use of cache for Jacobian (https://github.com/JuliaDiff/SparseDiffTools.jl/issues/135)
    function update_∇G(J, u, p, t)
        G_no_p(du, u) = G(du, u, p)
        forwarddiff_color_jacobian!(J, G_no_p, u, colorvec=colorsG, sparsity=jac0)
    end
    update_∇G(J, u, λ::Vector, t) = update_∇G(J, u, λ2p(λ), t)
    jacG = similar(jac0)
    ∇G = DiffEqArrayOperator(jacG, update_func=update_∇G)
    # ii) Create T as a BlockDiagDiffEqOperator and return jac = ∇G(x) - T
    T = BlockDiagDiffEqOperator{Float64}(Ts)
    jac = ∇G - T
    return ODEFunction{true}(f, jac=jac)
end
export AIBECSFunction

"""
    F, ∇ₓF = state_function_and_Jacobian(Ts, Gs)

Returns the state function `F` and its jacobian, `∇ₓF`.
"""
function F_and_∇ₓF(fun::ODEFunction, λ2p)
    ∇ₓF(u, p) = (update_coefficients!(fun.jac, u, p, nothing) ; convert(AbstractMatrix, fun.jac))
    ∇ₓF(u, λ::Vector) = ∇ₓF(u, λ2p(λ))
    F(du, u, p) = fun.f(du, u, p, nothing)
    F(u, p) = fun.f(u, p, nothing)
    F, ∇ₓF
end
F_and_∇ₓF(Ts, Gs, ::Type{P}) where {P <: APar} = F_and_∇ₓF(Ts, Gs, λ -> λ2p(P, λ))
F_and_∇ₓF(Ts, Gs, λ2p=nothing) = F_and_∇ₓF(AIBECSFunction(Ts, Gs, λ2p), λ2p)
export F_and_∇ₓF




"""
    localderivative(G, x, p)
    localderivative(Gᵢ, xs, i, p)
    localderivative(Gᵢ, dx, xs, i, p)

Returns the "local" derivative of `G` (or `Gᵢ`), i.e., equivalent to the vector

```
∇ₓG(x,p) * ones(size(x))
```

but using ForwardDiff's Jacobian instead.
"""
function localderivative(G, x, p) # for single tracer
    return ForwardDiff.derivative(λ -> G(x .+ λ, p), 0.0)
end
function localderivative(Gᵢ, xs, j, p) # for multiple tracers
    return ForwardDiff.derivative(λ -> Gᵢ(perturb_tracer(xs, j, λ)..., p), 0.0)
end
function localderivative(Gᵢ!, dx, xs, j, p) # if Gᵢ are in-place
    return ForwardDiff.derivative((dx, λ) -> Gᵢ!(dx, perturb_tracer(xs, j, λ)..., p), dx, 0.0)
end
function localderivative!(res, Gᵢ!, dx, xs, j, p) # if Gᵢ are in-place
    return ForwardDiff.derivative!(res, (dx, λ) -> Gᵢ!(dx, perturb_tracer(xs, j, λ)..., p), dx, 0.0)
end
perturb_tracer(xs, j, λ) = (xs[1:j - 1]..., xs[j] .+ λ, xs[j + 1:end]...)






# TODO update split functionality with DiffEqOperators
# This is important because it will eventually be used for
# preconditioning in Krylov methods
# Although I may need to revisit the API before that.
"""
    F, ∇ₓF = state_function_and_Jacobian(Ts, Gs, nb)

Returns the state function `F` and its jacobian, `∇ₓF`.
"""
function split_state_function_and_Jacobian(Ts::Tuple, Ls::Tuple, NLs::Tuple, nb)
    nt = length(Ts)
    tracers(x) = state_to_tracers(x, nb, nt)
    T(p) = blockdiag([Tⱼ(p) for Tⱼ in Ts]...) # Big T (linear part)
    NL(x, p) = reduce(vcat, NLⱼ(tracers(x)..., p) for NLⱼ in NLs) # nonlinear part
    L(x, p) = reduce(vcat, Lⱼ(tracers(x)..., p) for Lⱼ in Ls) # nonlinear part
    F(x, p) = NL(x, p) + L(x, p) - T(p) * x                     # full 𝐹(𝑥) = -T 𝑥 + 𝐺(𝑥)
    ∇ₓNL(x, p) = local_jacobian(NLs, x, p, nt, nb)     # Jacobian of nonlinear part
    ∇ₓL(p) = local_jacobian(Ls, zeros(nt * nb), p, nt, nb)     # Jacobian of nonlinear part
    ∇ₓF(x, p) = ∇ₓNL(x, p) + ∇ₓL(p) - T(p)       # full Jacobian ∇ₓ𝐹(𝑥) = -T + ∇ₓ𝐺(𝑥)
    return F, L, NL, ∇ₓF, ∇ₓL, ∇ₓNL, T
end
function split_state_function_and_Jacobian(T, L, NL, nb)
    F(x, p) = NL(x, p) + L(x, p) - T(p) * x                     # full 𝐹(𝑥)
    ∇ₓNL(x, p) = sparse(Diagonal(localderivative(NL, x, p))) # Jacobian of nonlinear part
    ∇ₓL(p) = sparse(Diagonal(localderivative(L, zeros(nb), p)))     # Jacobian of nonlinear part
    ∇ₓF(x, p) = ∇ₓNL(x, p) + ∇ₓL(p) - T(p)       # full Jacobian ∇ₓ𝐹(𝑥) = -T + ∇ₓ𝐺(𝑥)
    return F, L, NL, ∇ₓF, ∇ₓL, ∇ₓNL, T
end
export split_state_function_and_Jacobian

function local_jacobian(Gs, x, p, nt, nb)
    return reduce(vcat, local_jacobian_row(Gⱼ, x, p, nt, nb) for Gⱼ in Gs)
end
function local_jacobian2(Gs, x, p, nt, nb)
    return reduce(hcat, local_jacobian_col(Gs, x, p, nt, nb) for i in 1:nt)
end
function inplace_local_jacobian(Gs, x, p, nt, nb)
    return reduce(vcat, inplace_local_jacobian_row(Gⱼ!, x, p, nt, nb) for Gⱼ! in Gs)
end

function local_jacobian_row(Gᵢ, x, p, nt, nb)
    tracers(x) = state_to_tracers(x, nb, nt)
    return reduce(hcat, sparse(Diagonal(localderivative(Gᵢ, tracers(x), j, p))) for j in 1:nt)
end
function local_jacobian_col(Gs, x, p, nt, nb)
    tracers(x) = state_to_tracers(x, nb, nt)
    return reduce(vcat, sparse(Diagonal(localderivative2(Gᵢ, tracers(x), j, p))) for Gᵢ in 1:nt)
end
function inplace_local_jacobian_row(Gᵢ!, x, p, nt, nb)
    tracers(x) = state_to_tracers(x, nb, nt)
    dx = Vector{Float64}(undef, nb)
    return reduce(hcat, sparse(Diagonal(localderivative(Gᵢ!, dx, tracers(x), j, p))) for j in 1:nt)
end

#= ============================================
Generate 𝑓 and derivatives from user input
============================================ =#

function generate_f(ωs, μx, σ²x, v, ωp, ::Type{T}) where {T <: APar}
    nt, nb = length(ωs), length(v)
    tracers(x) = state_to_tracers(x, nb, nt)
    f(x, λorp) = ωp * mismatch(T, λorp) +
        sum([ωⱼ * mismatch(xⱼ, μⱼ, σⱼ², v) for (ωⱼ, xⱼ, μⱼ, σⱼ²) in zip(ωs, tracers(x), μx, σ²x)])
    return f
end
function generate_f(ωs, ωp, grd, obs, ::Type{T}; kwargs...) where {T <: APar}
    nt, nb = length(ωs), count(iswet(grd))
    tracers(x) = state_to_tracers(x, nb, nt)
    Ms = [interpolationmatrix(grd, obsⱼ) for obsⱼ in obs]
    cs = get(kwargs, :cs, (collect(identity for i in 1:nt)...,))
    f(x, λorp) = ωp * mismatch(T, λorp) +
        sum([ωⱼ * mismatch(xⱼ, grd, obsⱼ, M=Mⱼ, c=cⱼ) for (ωⱼ, xⱼ, obsⱼ, Mⱼ, cⱼ) in zip(ωs, tracers(x), obs, Ms, cs)])
    return f
end
function generate_f(ωs, ωp, grd, modify::Function, obs, ::Type{T}) where {T <: APar}
    nt, nb = length(ωs), count(iswet(grd))
    Ms = [interpolationmatrix(grd, obsⱼ) for obsⱼ in obs]
    iwets = [iswet(grd, obsⱼ) for obsⱼ in obs]
    function f(x, λorp)
        xs = unpack_tracers(x, grd)
        return ωp * mismatch(T, λorp) + sum([ωᵢ * indirectmismatch(xs, grd, modify, obs, i, Mᵢ, iwetᵢ) for (i, (ωᵢ, Mᵢ, iwetᵢ)) in enumerate(zip(ωs, Ms, iwets))])
    end
    return f
end




function generate_∇ₓf(ωs, μx, σ²x, v)
    nt, nb = length(ωs), length(v)
    tracers(x) = state_to_tracers(x, nb, nt)
    ∇ₓf(x) = reduce(hcat, ωⱼ * ∇mismatch(xⱼ, μⱼ, σⱼ², v) for (ωⱼ, xⱼ, μⱼ, σⱼ²) in zip(ωs, tracers(x), μx, σ²x))
    ∇ₓf(x, p) = ∇ₓf(x)
    return ∇ₓf
end
function generate_∇ₓf(ωs, grd, obs; kwargs...)
    nt, nb = length(ωs), count(iswet(grd))
    tracers(x) = state_to_tracers(x, nb, nt)
    Ms = [interpolationmatrix(grd, obsⱼ) for obsⱼ in obs]
    cs = get(kwargs, :cs, (collect(identity for i in 1:nt)...,))
    ∇ₓf(x) = reduce(hcat, ωⱼ * ∇mismatch(xⱼ, grd, obsⱼ, M=Mⱼ, c=cⱼ) for (ωⱼ, xⱼ, obsⱼ, Mⱼ, cⱼ) in zip(ωs, tracers(x), obs, Ms, cs))
    ∇ₓf(x, p) = ∇ₓf(x)
    return ∇ₓf
end
function generate_∇ₓf(ωs, grd, modify::Function, obs)
    nt, nb = length(ωs), count(iswet(grd))
    Ms = [interpolationmatrix(grd, obsⱼ) for obsⱼ in obs]
    iwets = [iswet(grd, obsⱼ) for obsⱼ in obs]
    function ∇ₓf(x)
        xs = unpack_tracers(x, grd)
        sum([ωᵢ * ∇indirectmismatch(unpack_tracers(x, grd), grd, modify, obs, i, Mᵢ, iwetᵢ) for (i, (ωᵢ, Mᵢ, iwetᵢ)) in enumerate(zip(ωs, Ms, iwets))])
    end
    ∇ₓf(x, p) = ∇ₓf(x)
    return ∇ₓf
end


function f_and_∇ₓf(ωs, μx, σ²x, v, ωp, ::Type{T}) where {T <: APar}
    generate_f(ωs, μx, σ²x, v, ωp, T), generate_∇ₓf(ωs, μx, σ²x, v)
end
function f_and_∇ₓf(ωs, ωp, grd, obs, ::Type{T}; kwargs...) where {T <: APar}
    generate_f(ωs, ωp, grd, obs, T; kwargs...), generate_∇ₓf(ωs, grd, obs; kwargs...)
end
function f_and_∇ₓf(ωs, ωp, grd, modify::Function, obs, ::Type{T}) where {T <: APar}
    generate_f(ωs, ωp, grd, modify, obs, T), generate_∇ₓf(ωs, grd, modify, obs)
end
export f_and_∇ₓf

"""
    mismatch(x, xobs, σ²xobs, v)

Volume-weighted mismatch of modelled tracer `x` against observed mean, `xobs`, given observed variance, `σ²xobs`, and volumes `v`.
"""
function mismatch(x, xobs, σ²xobs, v)
    δx = x - xobs
    W = Diagonal(v ./ σ²xobs)
    return 0.5 * transpose(δx) * W * δx / (transpose(xobs) * W * xobs)
end
mismatch(x, ::Missing, args...) = 0

"""
    ∇mismatch(x, xobs, σ²xobs, v)

Adjoint of the gradient of `mismatch(x, xobs, σ²xobs, v)`.
"""
function ∇mismatch(x, xobs, σ²xobs, v)
    δx = x - xobs
    W = Diagonal(v ./ σ²xobs)
    return transpose(W * δx) / (transpose(xobs) * W * xobs)
end
∇mismatch(x, ::Missing, args...) = transpose(zeros(length(x)))




## new functions for more generic obs packages
# TODO Add an optional function argument to transform the data before computingn the mismatch
# Example if for isotope tracers X where one ususally wants to minimize the mismatch in δ or ε.
function mismatch(x, grd::OceanGrid, obs; c=identity, W=I, M=interpolationmatrix(grd, obs.metadata), iwet=iswet(grd, obs))
    o = view(obs, iwet)
    δx = M * c(x) - o
    return 0.5 * transpose(δx) * W * δx / (transpose(o) * W * o)
end
mismatch(x, grd::OceanGrid, ::Missing; kwargs...) = 0
function ∇mismatch(x, grd::OceanGrid, obs; c=identity, W=I, M=interpolationmatrix(grd, obs.metadata), iwet=iswet(grd, obs))
    ∇c = Diagonal(ForwardDiff.derivative(λ -> c(x .+ λ), 0.0))
    o = view(obs, iwet)
    δx = M * c(x) - o
    return transpose(W * δx) * M * ∇c / (transpose(o) * W * o)
end
∇mismatch(x, grd::OceanGrid, ::Missing; kwargs...) = transpose(zeros(length(x)))

# In case the mismatch is not based on the tracer but on some function of it
function indirectmismatch(xs::Tuple, grd::OceanGrid, modify::Function, obs, i, M=interpolationmatrix(grd, obs[i].metadata), iwet=iswet(grd, obs[i]))
    x2 = modify(xs...)
    out = 0.0
    M = interpolationmatrix(grd, obs[i].metadata)
    iwet = iswet(grd, obs[i])
    o = view(obs[i], iwet)
    δx = M * x2[i] - o
    return 0.5 * transpose(δx) * δx / (transpose(o) * o)
end
function ∇indirectmismatch(xs::Tuple, grd::OceanGrid, modify::Function, obs, i, M=interpolationmatrix(grd, obs[i].metadata), iwet=iswet(grd, obs[i]))
    nt, nb = length(xs), length(iswet(grd))
    x2 = modify(xs...)
    o = view(obs[i], iwet)
    δx = M * x2[i] - o
    ∇modᵢ = ∇modify(modify, xs, i)
    return transpose(δx) * M * ∇modᵢ / (transpose(o) * o)
end
# TODO think of more efficient way to avoid recomputing ∇modify whole for each i
function ∇modify(modify, xs, i, j)
    return sparse(Diagonal(ForwardDiff.derivative(λ -> modify(perturb_tracer(xs, j, λ)...)[i], 0.0)))
end
∇modify(modify, xs, i) = reduce(hcat, ∇modify(modify, xs, i, j) for j in 1:length(xs))


#= ============================================
multi-tracer norm
============================================ =#

function volumeweighted_norm(nt, v)
    w = repeat(v, nt)
    return nrm(x) = transpose(x) * Diagonal(w) * x
end


#= ============================================
unpacking of multi-tracers
============================================ =#

state_to_tracers(x, nb, nt) = ntuple(i -> state_to_tracer(x, nb, nt, i), nt)
state_to_tracer(x, nb, nt, i) = view(x, tracer_indices(nb, nt, i))
function state_to_tracers(x, grd)
    nb = number_of_wet_boxes(grd)
    nt = Int(round(length(x) / nb))
    return state_to_tracers(x, nb, nt)
end
tracer_indices(nb, nt, i) = (i - 1) * nb + 1:i * nb
tracers_to_state(xs) = reduce(vcat, xs)
export state_to_tracers, state_to_tracer, tracers_to_state, tracer_indices
# Alias for better name
unpack_tracers = state_to_tracers
export unpack_tracers




