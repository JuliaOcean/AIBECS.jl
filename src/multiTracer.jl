



#=============================================
Generate 𝐹 and ∇ₓ𝐹 from user input
=============================================#

"""
    F, ∇ₓF = state_function_and_Jacobian(Ts, Gs, nb)

Returns the state function `F` and its jacobian, `∇ₓF`.
"""
function state_function_and_Jacobian(Ts::Tuple, Gs::Tuple, nb)
    nt = length(Ts)
    tracers(x) = state_to_tracers(x, nb, nt)
    tracer(x,i) = state_to_tracer(x, nb, nt, i)
    T(p) = blockdiag([Tⱼ(p) for Tⱼ in Ts]...) # Big T (linear part)
    function G(x,p)
        xs = tracers(x)
        return reduce(vcat, [Gⱼ(xs..., p) for Gⱼ in Gs]) # nonlinear part
    end
    F(x,p) = G(x,p) - reduce(vcat, [Tⱼ(p) * tracer(x,j) for (j,Tⱼ) in enumerate(Ts)])
    ∇ₓG(x,p) = local_jacobian(Gs, x, p, nt, nb)     # Jacobian of nonlinear part
    ∇ₓF(x,p) = ∇ₓG(x,p) - T(p)       # full Jacobian ∇ₓ𝐹(𝑥) = -T + ∇ₓ𝐺(𝑥)
    return F, ∇ₓF
end
function state_function_and_Jacobian(T, G)
    F(x,p) = G(x,p) - T(p) * x                     # full 𝐹(𝑥) = -T 𝑥 + 𝐺(𝑥)
    ∇ₓG(x,p) = sparse(Diagonal(localderivative(G, x, p))) # Jacobian of nonlinear part
    ∇ₓF(x,p) = ∇ₓG(x,p) - T(p)       # full Jacobian ∇ₓ𝐹(𝑥) = -T + ∇ₓ𝐺(𝑥)
    return F, ∇ₓF
end

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
    return ForwardDiff.derivative(λ -> Gᵢ(perturb_tracer(xs,j,λ)..., p), 0.0)
end
function localderivative(Gᵢ!, dx, xs, j, p) # if Gᵢ are in-place
    return ForwardDiff.derivative((dx,λ) -> Gᵢ!(dx, perturb_tracer(xs,j,λ)..., p), dx, 0.0)
end
perturb_tracer(xs, j, λ) = (xs[1:j-1]..., xs[j] .+ λ, xs[j+1:end]...)

function inplace_state_function_and_Jacobian(Ts::Tuple, Gs::Tuple, nb)
    nt = length(Ts)
    tracers(x) = state_to_tracers(x, nb, nt)
    tracer(x,i) = state_to_tracer(x, nb, nt, i)
    T(p) = blockdiag([Tⱼ(p) for Tⱼ in Ts]...) # Big T (linear part)
    #F(x,p) = G(x,p) - T(p) * x                     # full 𝐹(𝑥) = -T 𝑥 + 𝐺(𝑥)
    function F!(dx,x,p)
        xs = tracers(x)
        for (j, (Tⱼ, Gⱼ!)) in enumerate(zip(Ts, Gs))
            ij = tracer_indices(nb,nt,j)
            @views dx[ij] .= Gⱼ!(dx[ij], xs..., p)
            @views dx[ij] .-= Tⱼ(p) * x[ij]
        end
        return dx
    end
    ∇ₓG(x,p) = inplace_local_jacobian(Gs, x, p, nt, nb)     # Jacobian of nonlinear part
    ∇ₓF(x,p) = ∇ₓG(x,p) - T(p)       # full Jacobian ∇ₓ𝐹(𝑥) = -T + ∇ₓ𝐺(𝑥)
    return F!, ∇ₓF
end


export state_function_and_Jacobian, inplace_state_function_and_Jacobian

"""
    F, ∇ₓF = state_function_and_Jacobian(Ts, Gs, nb)

Returns the state function `F` and its jacobian, `∇ₓF`.
"""
function split_state_function_and_Jacobian(Ts::Tuple, Ls::Tuple, NLs::Tuple, nb)
    nt = length(Ts)
    tracers(x) = state_to_tracers(x, nb, nt)
    T(p) = blockdiag([Tⱼ(p) for Tⱼ in Ts]...) # Big T (linear part)
    NL(x,p) = reduce(vcat, [NLⱼ(tracers(x)..., p) for NLⱼ in NLs]) # nonlinear part
    L(x,p) = reduce(vcat, [Lⱼ(tracers(x)..., p) for Lⱼ in Ls]) # nonlinear part
    F(x,p) = NL(x,p) + L(x,p) - T(p) * x                     # full 𝐹(𝑥) = -T 𝑥 + 𝐺(𝑥)
    ∇ₓNL(x,p) = local_jacobian(NLs, x, p, nt, nb)     # Jacobian of nonlinear part
    ∇ₓL(p) = local_jacobian(Ls, zeros(nt*nb), p, nt, nb)     # Jacobian of nonlinear part
    ∇ₓF(x,p) = ∇ₓNL(x,p) + ∇ₓL(p) - T(p)       # full Jacobian ∇ₓ𝐹(𝑥) = -T + ∇ₓ𝐺(𝑥)
    return F, L, NL, ∇ₓF, ∇ₓL, ∇ₓNL, T
end
function split_state_function_and_Jacobian(T, L, NL, nb)
    F(x,p) = NL(x,p) + L(x,p) - T(p) * x                     # full 𝐹(𝑥)
    ∇ₓNL(x,p) = sparse(Diagonal(localderivative(NL, x, p))) # Jacobian of nonlinear part
    ∇ₓL(p) = sparse(Diagonal(localderivative(L, zeros(nb), p)))     # Jacobian of nonlinear part
    ∇ₓF(x,p) = ∇ₓNL(x,p) + ∇ₓL(p) - T(p)       # full Jacobian ∇ₓ𝐹(𝑥) = -T + ∇ₓ𝐺(𝑥)
    return F, L, NL, ∇ₓF, ∇ₓL, ∇ₓNL, T
end
export split_state_function_and_Jacobian

function local_jacobian(Gs, x, p, nt, nb)
    return reduce(vcat, [local_jacobian_row(Gⱼ, x, p, nt, nb) for Gⱼ in Gs])
end
function inplace_local_jacobian(Gs, x, p, nt, nb)
    return reduce(vcat, [inplace_local_jacobian_row(Gⱼ!, x, p, nt, nb) for Gⱼ! in Gs])
end

function local_jacobian_row(Gᵢ, x, p, nt, nb)
    tracers(x) = state_to_tracers(x, nb, nt)
    return reduce(hcat, [sparse(Diagonal(localderivative(Gᵢ, tracers(x), j, p))) for j in 1:nt])
end
function inplace_local_jacobian_row(Gᵢ!, x, p, nt, nb)
    tracers(x) = state_to_tracers(x, nb, nt)
    dx = Vector{Float64}(undef,nb)
    return reduce(hcat, [sparse(Diagonal(localderivative(Gᵢ!, dx, tracers(x), j, p))) for j in 1:nt])
end

#=============================================
Generate 𝑓 and derivatives from user input
=============================================#

function generate_objective(ωs, μx, σ²x, v, ωp)
    nt, nb = length(ωs), length(v)
    tracers(x) = state_to_tracers(x, nb, nt)
    f(x, p) = ωp * mismatch(p) +
        sum([ωⱼ * mismatch(xⱼ, μⱼ, σⱼ², v) for (ωⱼ, xⱼ, μⱼ, σⱼ²) in zip(ωs, tracers(x), μx, σ²x)])
    return f
end
function generate_objective(ωs, ωp, grd, obs; kwargs...)
    nt, nb = length(ωs), count(iswet(grd))
    tracers(x) = state_to_tracers(x, nb, nt)
    Ms = [interpolationmatrix(grd, obsⱼ) for obsⱼ in obs]
    cs = get(kwargs, :cs, (collect(identity for i in 1:nt)...,))
    f(x, p) = ωp * mismatch(p) +
        sum([ωⱼ * mismatch(xⱼ, grd, obsⱼ, M=Mⱼ, c=cⱼ) for (ωⱼ, xⱼ, obsⱼ, Mⱼ, cⱼ) in zip(ωs, tracers(x), obs, Ms, cs)])
    return f
end
function generate_objective(ωs, ωp, grd, modify::Function, obs)
    nt, nb = length(ωs), count(iswet(grd))
    Ms = [interpolationmatrix(grd, obsⱼ) for obsⱼ in obs]
    iwets = [iswet(grd, obsⱼ) for obsⱼ in obs]
    function f(x, p)
        xs = unpack_tracers(x, grd)
        return ωp * mismatch(p) + sum([ωᵢ * indirectmismatch(xs, grd, modify, obs, i, Mᵢ, iwetᵢ) for (i, (ωᵢ, Mᵢ, iwetᵢ)) in enumerate(zip(ωs, Ms, iwets))])
    end
    return f
end




function generate_∇ₓobjective(ωs, μx, σ²x, v)
    nt, nb = length(ωs), length(v)
    tracers(x) = state_to_tracers(x, nb, nt)
    ∇ₓf(x, p) = reduce(hcat, [ωⱼ * ∇mismatch(xⱼ, μⱼ, σⱼ², v) for (ωⱼ, xⱼ, μⱼ, σⱼ²) in zip(ωs, tracers(x), μx, σ²x)])
    return ∇ₓf
end
function generate_∇ₓobjective(ωs, grd, obs; kwargs...)
    nt, nb = length(ωs), count(iswet(grd))
    tracers(x) = state_to_tracers(x, nb, nt)
    Ms = [interpolationmatrix(grd, obsⱼ) for obsⱼ in obs]
    cs = get(kwargs, :cs, (collect(identity for i in 1:nt)...,))
    ∇ₓf(x, p) = reduce(hcat, [ωⱼ * ∇mismatch(xⱼ, grd, obsⱼ, M=Mⱼ, c=cⱼ) for (ωⱼ, xⱼ, obsⱼ, Mⱼ, cⱼ) in zip(ωs, tracers(x), obs, Ms, cs)])
    return ∇ₓf
end
function generate_∇ₓobjective(ωs, grd, modify::Function, obs)
    nt, nb = length(ωs), count(iswet(grd))
    Ms = [interpolationmatrix(grd, obsⱼ) for obsⱼ in obs]
    iwets = [iswet(grd, obsⱼ) for obsⱼ in obs]
    function ∇ₓf(x, p)
        xs = unpack_tracers(x, grd)
        return sum([ωᵢ * ∇indirectmismatch(unpack_tracers(x, grd), grd, modify, obs, i, Mᵢ, iwetᵢ) for (i, (ωᵢ, Mᵢ, iwetᵢ)) in enumerate(zip(ωs, Ms, iwets))])
    end
    return ∇ₓf
end


function generate_∇ₚobjective(ωp)
    ∇ₚf(x, p) = ωp * ∇mismatch(p)
    return ∇ₚf
end


generate_objective_and_derivatives(ωs, μx, σ²x, v, ωp) =
    generate_objective(ωs, μx, σ²x, v, ωp),
    generate_∇ₓobjective(ωs, μx, σ²x, v),
    generate_∇ₚobjective(ωp)
generate_objective_and_derivatives(ωs, ωp, grd, obs; kwargs...) =
    generate_objective(ωs, ωp, grd, obs; kwargs...),
    generate_∇ₓobjective(ωs, grd, obs; kwargs...),
    generate_∇ₚobjective(ωp)
function generate_objective_and_derivatives(ωs, ωp, grd, modify::Function, obs)
    return generate_objective(ωs, ωp, grd, modify, obs),
             generate_∇ₓobjective(ωs, grd, modify, obs),
             generate_∇ₚobjective(ωp)
end

export generate_objective, generate_∇ₓobjective, generate_∇ₚobjective
export generate_objective_and_derivatives



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
function indirectmismatch(xs::Tuple, grd::OceanGrid, modify::Function, obs, i, M = interpolationmatrix(grd, obs[i].metadata), iwet = iswet(grd, obs[i]))
    x2 = modify(xs...)
    out = 0.0
    M = interpolationmatrix(grd, obs[i].metadata)
    iwet = iswet(grd, obs[i])
    o = view(obs[i], iwet)
    δx = M * x2[i] - o
    return 0.5 * transpose(δx) * δx / (transpose(o) * o)
end
function ∇indirectmismatch(xs::Tuple, grd::OceanGrid, modify::Function, obs, i, M = interpolationmatrix(grd, obs[i].metadata), iwet = iswet(grd, obs[i]))
    nt, nb = length(xs), length(iswet(grd))
    x2 = modify(xs...)
    o = view(obs[i], iwet)
    δx = M * x2[i] - o
    ∇modᵢ = ∇modify(modify, xs, i)
    return transpose(δx) * M * ∇modᵢ / (transpose(o) * o)
end
# TODO think of more efficient way to avoid recomputing ∇modify whole for each i
function ∇modify(modify, xs, i, j)
    return sparse(Diagonal(ForwardDiff.derivative(λ -> modify(perturb_tracer(xs,j,λ)...)[i], 0.0)))
end
∇modify(modify, xs, i) = reduce(hcat, [∇modify(modify, xs, i, j) for j in 1:length(xs)])


#=============================================
multi-tracer norm
=============================================#

function volumeweighted_norm(nt, v)
    w = repeat(v, nt)
    return nrm(x) = transpose(x) * Diagonal(w) * x
end


#=============================================
unpacking of multi-tracers
=============================================#

state_to_tracers(x, nb, nt) = ntuple(i -> state_to_tracer(x, nb, nt, i), nt)
state_to_tracer(x, nb, nt, i) = x[tracer_indices(nb, nt, i)]
function state_to_tracers(x, grd)
    nb = number_of_wet_boxes(grd)
    nt = Int(round(length(x) / nb))
    return state_to_tracers(x, nb, nt)
end
tracer_indices(nb, nt, i) = (i-1)*nb+1 : i*nb
tracers_to_state(xs) = reduce(vcat, xs)
export state_to_tracers, state_to_tracer, tracers_to_state, tracer_indices
# Alias for better name
unpack_tracers = state_to_tracers
export unpack_tracers




