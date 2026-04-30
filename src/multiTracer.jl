



#= ============================================
Generate 𝐹 and ∇ₓ𝐹 from user input
============================================ =#

# TODO replace this with DiffEqOperators when possible
"""
    LinearOperators(ops...)

Container for a tuple of linear operators that acts as their sum.

Lets AIBECS represent a transport operator built from several pieces (e.g. a
sparse advection matrix plus a matrix-free diffusion operator) without
materialising the sum eagerly. Supports `*`, `\\`, `+` with sparse matrices,
and in-place `mul!`. `factorize`, `*`, and `\\` materialise the sum on demand.

```julia
A = LinearOperators(T1, T2)   # behaves like T1 + T2
A * x                          # equivalent to (T1 + T2) * x
A \\ b                         # equivalent to (T1 + T2) \\ b
```
"""
struct LinearOperators{T<:Tuple}
    ops::T
end
LinearOperators(ops...) = LinearOperators((ops...,))
SparseArrays.blockdiag(As::LinearOperators...) = blockdiag([sum(A.ops) for A in As]...)
function LinearAlgebra.mul!(du, A::LinearOperators, u, α, β)
    for (i, op) in enumerate(A.ops)
        (i==1) ? mul!(du, op, u, α, β) : mul!(du, op, u, α, 1)
    end
    du
end
LinearAlgebra.factorize(A::LinearOperators) = factorize(sum(A.ops))
Base.:\(A::LinearOperators, u::AbstractArray) = factorize(A) \ u
Base.:+(A::LinearOperators, B::AbstractSparseArray) = LinearOperators(A.ops..., B)
Base.:+(B::AbstractSparseArray, A::LinearOperators) = A + B
function Base.:*(A::LinearOperators, u::AbstractArray) # has to be same type!
    du = similar(u)
    mul!(du, A, u, 1, 0)
end
export LinearOperators

"""
    AIBECSFunction(T, G, nb)
    AIBECSFunction(Ts, Gs, nb, [nt, Tidx])
    AIBECSFunction(Ts, Gs, nb, ::Type{P}) where {P <: AbstractParameters}

Build a `SciMLBase.ODEFunction` for the AIBECS state equation

```math
\\partial_t \\boldsymbol{x} = \\boldsymbol{G}(\\boldsymbol{x}, \\boldsymbol{p}) - \\mathbf{T}(\\boldsymbol{p}) \\, \\boldsymbol{x}
```

together with its Jacobian ``\\nabla_x F = \\nabla_x \\boldsymbol{G} - \\mathbf{T}``,
ready to plug into `SteadyStateProblem` or any SciML solver.

`T` is the linear transport operator: a sparse matrix (single tracer, constant
in `p`) or a function `p -> T(p)` (depends on parameters). `G` is the local
sources-and-sinks function `G(x, p)`. For multi-tracer models pass tuples `Ts`
and `Gs`; the optional `Tidx` maps each of the `nt` tracers to one of the
operators in `Ts` (defaults to one operator per tracer). When a parameter type
`P <: AbstractParameters` is supplied, the returned function also accepts a
flat `Vector` parameter (converted internally via [`λ2p`](@ref)).

```julia
T(p) = transportoperator(grd, OCIM2.load()[2])
G(x, p) = -x ./ p.τ
fun = AIBECSFunction(T, G, count(iswet(grd)))
prob = SteadyStateProblem(fun, x₀, p)
```
"""
AIBECSFunction(T::AbstractSparseArray, G::Function) = AIBECSFunction(p -> T, G, T.n)
AIBECSFunction(T::Function, G::Function, nb::Int) = AIBECSFunction(T, (G,), nb)
function AIBECSFunction(T::Function, Gs::Tuple, nb::Int)
    nt = length(Gs) # if a single T is given, then nt is given by the number of Gs
    Tidx = ones(Int64, nt) # and all the tracers share the same T
    AIBECSFunction((T,), Gs, nb, nt, Tidx)
end
function AIBECSFunction(Ts::Tuple, Gs::Tuple, nb::Int, nt::Int=length(Ts), Tidx::AbstractVector=1:nt)
    tracers(u) = state_to_tracers(u, nb, nt)
    tracer(u, i) = state_to_tracer(u, nb, nt, i)
    G(u, p) = reduce(vcat, Gⱼ(tracers(u)..., p) for Gⱼ in Gs)
    function f(u, p, t=0)
        du = copy(G(u, p))
        for jT in eachindex(Ts)
            op = Ts[jT](p)
            for j in findall(Tidx .== jT)
                mul!(tracer(du, j), op, tracer(u, j), -1, 1)
            end
        end
        du
    end
    # Jacobian
    ∇ₓG(u, p) = local_jacobian(Gs, u, p, nt, nb)
    function T(p)
        uniqueTs = [Tⱼ(p) for Tⱼ in Ts]
        blockdiag([uniqueTs[Tidx[j]] for j in 1:nt]...)
    end
    jac(u, p, t=0) = ∇ₓG(u, p) - T(p)
    return ODEFunction(f, jac=jac)
end
# AIBECSFunction calls itself to allow for both λ::Vector and p::APar
function AIBECSFunction(Ts, Gs, nb::Int, ::Type{P}) where {P <: APar}
    AIBECSFunction(AIBECSFunction(Ts, Gs, nb), P)
end
function AIBECSFunction(fun::ODEFunction, ::Type{P}) where {P <: APar}
    jac(u, p::P, t=0) = fun.jac(u, p, t)
    jac(u, λ::Vector, t=0) = fun.jac(u, λ2p(P, λ), t)
    f(u, p::P, t=0) = fun.f(u, p, t)
    f(u, λ::Vector, t=0) = fun.f(u, λ2p(P, λ), t)
    return ODEFunction{false}(f, jac=jac)
end

export AIBECSFunction

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
perturb_tracer(xs, j, λ) = (xs[1:j - 1]..., xs[j] .+ λ, xs[j + 1:end]...)



"""
    F, L, NL, ∇ₓF, ∇ₓL, ∇ₓNL, T = split_state_function_and_Jacobian(Ts, Ls, NLs, nb)

Returns the state function `F` and its jacobian, `∇ₓF`, as well as a collection of split operators. This is experimental. Use at your own risk!
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

function local_jacobian_row(Gᵢ, x, p, nt, nb)
    tracers(x) = state_to_tracers(x, nb, nt)
    return reduce(hcat, sparse(Diagonal(localderivative(Gᵢ, tracers(x), j, p))) for j in 1:nt)
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
function mismatch(x, grd::OceanGrid, obs; c=identity, W=I, M=interpolationmatrix(grd, obs), iwet=iswet(grd, obs))
    o = view(obs, iwet)
    δx = M * c(x) - o
    return 0.5 * transpose(δx) * W * δx / (transpose(o) * W * o)
end
mismatch(x, grd::OceanGrid, ::Missing; kwargs...) = 0
function ∇mismatch(x, grd::OceanGrid, obs; c=identity, W=I, M=interpolationmatrix(grd, obs), iwet=iswet(grd, obs))
    ∇c = Diagonal(ForwardDiff.derivative(λ -> c(x .+ λ), 0.0))
    o = view(obs, iwet)
    δx = M * c(x) - o
    return transpose(W * δx) * M * ∇c / (transpose(o) * W * o)
end
∇mismatch(x, grd::OceanGrid, ::Missing; kwargs...) = transpose(zeros(length(x)))

# In case the mismatch is not based on the tracer but on some function of it
# TODO Add option for correlation matrix for Bayesian inferenece
function indirectmismatch(xs::Tuple, grd::OceanGrid, modify::Function, obs, i, M=interpolationmatrix(grd, obs[i]), iwet=iswet(grd, obs[i]))
    x2 = modify(xs...)
    out = 0.0
    o = obs[i][iwet, :value]
    δx = M * x2[i] - o
    return 0.5 * transpose(δx) * δx / (transpose(o) * o)
end
function ∇indirectmismatch(xs::Tuple, grd::OceanGrid, modify::Function, obs, i, M=interpolationmatrix(grd, obs[i]), iwet=iswet(grd, obs[i]))
    x2 = modify(xs...)
    o = obs[i][iwet, :value]
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

"""
    state_to_tracers(x, nb, nt)
    state_to_tracers(x, grd)

Split the flat state vector `x` into a tuple of `nt` tracer views of length `nb`.

The returned objects are `view`s into `x`, so mutating them mutates `x`. The
2-argument form infers `nb` from `grd` and `nt` from `length(x) / nb`.

```julia
xs = state_to_tracers(x, grd)
DIP, DOP = xs                    # for a 2-tracer model
```
"""
state_to_tracers(x, nb, nt) = ntuple(i -> state_to_tracer(x, nb, nt, i), nt)

"""
    state_to_tracer(x, nb, nt, i)

View into the `i`th tracer of the flat state vector `x` (length `nb*nt`).
Equivalent to `state_to_tracers(x, nb, nt)[i]`.
"""
state_to_tracer(x, nb, nt, i) = view(x, tracer_indices(nb, nt, i))
function state_to_tracers(x, grd)
    nb = number_of_wet_boxes(grd)
    nt = Int(round(length(x) / nb))
    return state_to_tracers(x, nb, nt)
end

"""
    tracer_indices(nb, nt, i)

Range of indices in the flat state vector that correspond to the `i`th tracer
of an `nt`-tracer, `nb`-box model.
"""
tracer_indices(nb, nt, i) = (i - 1) * nb + 1:i * nb

"""
    tracers_to_state(xs)

Concatenate a tuple/array of tracer vectors `xs` into a flat state vector,
the inverse of [`state_to_tracers`](@ref).
"""
tracers_to_state(xs) = reduce(vcat, xs)
export state_to_tracers, state_to_tracer, tracers_to_state, tracer_indices

"""
    unpack_tracers(x, grd)
    unpack_tracers(x, nb, nt)

Alias for [`state_to_tracers`](@ref) with a more user-friendly name.
"""
unpack_tracers(args...) = state_to_tracers(args...)
export unpack_tracers




