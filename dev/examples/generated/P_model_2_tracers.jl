using AIBECS

grd, T_OCIM = OCIM0.load() ;

T_DIP(p) = T_OCIM

T_POP(p) = buildPFD(grd, settling_velocity = w(p))

w(p) = p.w₀ .+ p.w′ * z

iwet = findall(vec(grd.wet3D))
z = ustrip.(grd.depth_3D[iwet])

T_all = (T_DIP, T_POP) ;

geores(x, p) = (p.DIPgeo .- x) / p.τgeo

relu(x) = (x .≥ 0) .* x

function uptake(DIP, p)
    τ, k, z₀ = p.τ, p.k, p.z₀
    DIP⁺ = relu(DIP)
    return 1/τ * DIP⁺.^2 ./ (DIP⁺ .+ k) .* (z .≤ z₀)
end

remineralization(POP, p) = p.κ * POP

function sms_DIP(DIP, POP, p)
    return -uptake(DIP, p) + remineralization(POP, p) + geores(DIP, p)
end
function sms_POP(DIP, POP, p)
    return  uptake(DIP, p) - remineralization(POP, p)
end

sms_all = (sms_DIP, sms_POP) # bundles all the source-sink functions in a tuple

function G_DIP!(dDIP, DIP, POP, p)
    τ, k, z₀, κ, DIPgeo, τgeo = p.τ, p.k, p.z₀, p.κ, p.DIPgeo, p.τgeo
    dDIP .= @. -(DIP ≥ 0) / τ * DIP^2 / (DIP + k) * (z ≤ z₀) + κ * POP + (DIPgeo - DIP) / τgeo
    return dDIP
end

function G_POP!(dPOP, DIP, POP, p)
    τ, k, z₀, κ = p.τ, p.k, p.z₀, p.κ
    dPOP .= @. (DIP ≥ 0) / τ * DIP^2 / (DIP + k) * (z ≤ z₀) - κ * POP
    return dPOP
end

Gs = (G_DIP!, G_POP!)

t = empty_parameter_table()    # empty table of parameters

add_parameter!(t, :DIPgeo, 2.12u"mmol/m^3",
               optimizable = true,
               variance_obs = ustrip(upreferred(0.1 * 2.17u"mmol/m^3"))^2,
               description = "Mean PO₄ concentration")
add_parameter!(t, :τgeo, 1.0u"Myr",
               description = "Geological restoring timescale")
add_parameter!(t, :k, 6.62u"μmol/m^3",
               optimizable = true,
               description = "Half-saturation constant (Michaelis-Menten)")
add_parameter!(t, :z₀, 80.0u"m",
               description = "Depth of the euphotic layer base")
add_parameter!(t, :w₀, 0.64u"m/d",
               optimizable = true,
               description = "Sinking velocity at surface")
add_parameter!(t, :w′, 0.13u"1/d",
               optimizable = true,
               description = "Vertical gradient of sinking velocity")
add_parameter!(t, :κ, 0.19u"1/d",
               optimizable = true,
               description = "Remineralization rate constant (POP to DIP)")
add_parameter!(t, :τ, 236.52u"d",
               optimizable = true,
               description = "Maximum uptake rate timescale")

initialize_Parameters_type(t, "Pcycle_Parameters")   # Generate the parameter type

p = Pcycle_Parameters()

nb = length(iwet)
F!, ∇ₓF = inplace_state_function_and_Jacobian(T_all, Gs, nb)
F(x::Vector{Tx}, p::Pcycle_Parameters{Tp}) where {Tx,Tp} = F!(Vector{promote_type(Tx,Tp)}(undef,length(x)),x,p)

x = p.DIPgeo * ones(2nb) # initial iterate

prob = SteadyStateProblem(F, ∇ₓF, x, p)

s = solve(prob, CTKAlg()).u

##md # !!! warn

using Plots

DIP, POP = state_to_tracers(s, grd)

horizontalslice(DIP, grd, 2000; color=:viridis)

using WorldOceanAtlasTools
μDIPobs3D, σ²DIPobs3D = WorldOceanAtlasTools.fit_to_grid(grd, 2018, "phosphate", "annual", "1°", "an")
μDIPobs, σ²DIPobs = μDIPobs3D[iwet], σ²DIPobs3D[iwet]
μx = (μDIPobs, missing)
σ²x = (σ²DIPobs, missing)

ωs = [1.0, 0.0] # the weight for the mismatch (weight of POP = 0)
ωp = 1e-4       # the weight for the parameters prior estimates

v = ustrip.(vector_of_volumes(grd))
f = generate_objective(ωs, μx, σ²x, v, ωp, mean_obs(p), variance_obs(p))
∇ₓf = generate_∇ₓobjective(ωs, μx, σ²x, v, ωp, mean_obs(p), variance_obs(p))
∇ₚf = generate_∇ₚobjective(ωs, μx, σ²x, v, ωp, mean_obs(p), variance_obs(p))

using F1Method

mem = F1Method.initialize_mem(s, p)

objective(p) = F1Method.objective(f, F, ∇ₓF, mem, p, CTKAlg())
gradient(p) = F1Method.gradient(f, F, ∇ₓf, ∇ₓF, mem, p, CTKAlg())
hessian(p) = F1Method.hessian(f, F, ∇ₓf, ∇ₓF, mem, p, CTKAlg())

λ2p(λ) = AIBECS.opt_para(exp.(λ))

∇λ2p(λ) = exp.(λ)
∇²λ2p(λ) = exp.(λ)

p2λ(p) = log.(optvec(p))
λ = p2λ(p)

obj(λ) = objective(λ2p(λ))
using LinearAlgebra
grad(λ) = gradient(λ2p(λ)) * Diagonal(∇λ2p(λ))
function hess(λ)
    ∇p = Diagonal(∇λ2p(λ)) # for variable change
    ∇²p = Diagonal(∇²λ2p(λ)) # for variable change
    G = vec(gradient(λ2p(λ)))
    H = hessian(λ2p(λ))
    return ∇p * H * ∇p + Diagonal(G) * ∇²p
end

using Optim

m = length(p)
grad(s, λ) = s[1:m] .= vec(grad(λ))
hess(s, λ) = s[1:m,1:m] .= hess(λ)

opt = Optim.Options(store_trace = false, show_trace = true, extended_trace = false)

p

results = optimize(obj, grad, hess, λ, NewtonTrustRegion(), opt)

p_optimized = λ2p(results.minimizer)

prob_optimized = SteadyStateProblem(F, ∇ₓF, s, p_optimized)
s_optimized = solve(prob_optimized, CTKAlg()).u

optimized_DIP, POPnew = state_to_tracers(s_optimized, grd)
δDIP_before = 100(DIP - μDIPobs) ./ μDIPobs
δDIP_after = 100(optimized_DIP - μDIPobs) ./ μDIPobs

horizontalslice(δDIP_before, grd, 2000; levels=-20:2:20, color=:balance)

horizontalslice(δDIP_after, grd, 2000; levels=-20:2:20, color=:balance)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

