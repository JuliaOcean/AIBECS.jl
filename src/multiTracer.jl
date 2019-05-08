
#=============================================
Generate 𝐹 and ∇ₓ𝐹 from user input
=============================================#

# Create F and ∇ₓF automatically from Ts and Gs only
function state_function_and_Jacobian(Ts, Gs, nt, nb)
    tracers(v) = [v[j:j+nb-1] for j in 1:nb:nb*nt]
    T(p) = blockdiag([Tⱼ(p) for Tⱼ in Ts]...) # Big T (linear part)
    G(x, p) = reduce(vcat, [Gⱼ(tracers(x)..., p) for Gⱼ in Gs]) # nonlinear part
    F(x, p) = -T(p) * x + G(x, p)                     # full 𝐹(𝑥) = T 𝑥 + 𝐺(𝑥)
    ∇ₓG(x, p) = local_jacobian(Gs, x, p, nt, nb)     # Jacobian of nonlinear part
    ∇ₓF(x, p) = -T(p) + ∇ₓG(x, p)          # full Jacobian ∇ₓ𝐹(𝑥) = T + ∇ₓ𝐺(𝑥)
    return F, ∇ₓF
end
export state_function_and_Jacobian

function local_jacobian(Gs, x, p, nt, nb)
    return reduce(vcat, [local_jacobian_row(Gⱼ, x, p, nt, nb) for Gⱼ in Gs])
end

𝔇(x) = DualNumbers.dualpart.(x)      # dual part

function local_jacobian_row(Gⱼ, x, p, nt, nb)
    e(j) = kron([j == k for k in 1:nt], trues(nb))
    tracers(v) = [v[j:j+nb-1] for j in 1:nb:nb*nt]
    return reduce(hcat, [spdiagm(0 => 𝔇(Gⱼ(tracers(x + ε * e(j))..., p))) for j in 1:nt])
end

#=============================================
Generate 𝑓 and ∇ₓ𝑓 from user input
=============================================#

function mismatch_function_and_Jacobian(ωs, μx, σ²x, v, ωp, μp, σ²p)
    tracers(x) = [x[j:j+nb-1] for j in 1:nb:nb*nt]
    f(x, p) = ωp * mismatch(p, μp, σ²p) +
        sum([ωⱼ * mismatch(xⱼ, μⱼ, σⱼ², v) for (ωⱼ, xⱼ, μⱼ, σⱼ²) in zip(ωs, tracers(x), μx, σ²x)])
    ∇ₓf(x, p) = reduce(hcat, [ωⱼ * ∇mismatch(xⱼ, μⱼ, σⱼ², v) for (ωⱼ, xⱼ, μⱼ, σⱼ²) in zip(ωs, tracers(x), μx, σ²x)])
    return f, ∇ₓf
end
export mismatch_function_and_Jacobian

"""
    mismatch(x, xobs, σ²xobs, v)

Volume-weighted mismatch of modelled tracer `x` against observed mean, `xobs`, given observed variance, `σ²xobs`, and volumes `v`.
"""
function mismatch(x, xobs, σ²xobs, v)
    δx = x - xobs
    W = Diagonal(v ./ σ²xobs)
    return 0.5 * δx' * W * δx / (xobs' * W * xobs)
end

"""
    ∇mismatch(x, xobs, σ²xobs, v)

Adjoint of the gradient of `mismatch(x, xobs, σ²xobs, v)`.
"""
function ∇mismatch(x, xobs, σ²xobs, v)
    δx = x - xobs
    W = Diagonal(v ./ σ²xobs)
    return (W * δx)' / (xobs' * W * xobs)
end

# TODO
# Talk about it with FP
"""
    mismatch(p, logpobs, σ²logpobs)

    Mismatch of the log of model parameters `p` against observed (log) mean, `logpobs`, given observed (log) variance, `σ²logpobs`.
"""
function mismatch(p, logpobs, σ²logpobs)
    println("parameter mismatch to be checked!")
    δλ = log.(p) - logpobs
    W = Diagonal(1 ./ σ²logpobs)
    return 0.5 * δλ' * W * δλ
end

#=============================================
Generate multi-tracer norm
=============================================#

function volumeweighted_norm(nt, v)
    w = repeat(v, nt)
    return nrm(x) = transpose(x) * Diagonal(w) * x
end




