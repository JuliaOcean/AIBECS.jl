using AIBECS

wet3D, grd, T_OCIM = OCIM0.load() ;

T_DIP(p) = T_OCIM

T_POP(p) = buildPFD(grd, wet3D, sinking_speed = w(p))

w(p) = p.w₀ .+ p.w′ * z

iwet = findall(wet3D)
z = ustrip.(grd.depth_3D[iwet])

T_all = (T_DIP, T_POP) ;

geores(x, p) = (p.xgeo .- x) / p.τg

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

t = empty_parameter_table()    # empty table of parameters

add_parameter!(t, :xgeo, 2.12u"mmol/m^3",
               optimizable = true,
               variance_obs = ustrip(upreferred(0.1 * 2.17u"mmol/m^3"))^2,
               description = "Mean PO₄ concentration")
add_parameter!(t, :τg, 1.0u"Myr",
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
F, ∇ₓF = state_function_and_Jacobian(T_all, sms_all, nb)

x = p.xgeo * ones(2nb) # initial iterate

prob = SteadyStateProblem(F, ∇ₓF, x, p)

s = solve(prob, CTKAlg()).u

DIP, POP = state_to_tracers(s, nb, 2)

iz = findfirst(grd.depth .> 2000u"m")
iz, grd.depth[iz]

DIP_3D = rearrange_into_3Darray(DIP, wet3D)
DIP_2D = DIP_3D[:,:,iz] * ustrip(1.0u"mol/m^3" |> u"mmol/m^3")
lat, lon = ustrip.(grd.lat), ustrip.(grd.lon)

ENV["MPLBACKEND"]="qt5agg"
using PyPlot, PyCall
clf()
ccrs = pyimport("cartopy.crs")
ax = subplot(projection = ccrs.EqualEarth(central_longitude=-155.0))
ax.coastlines()

lon_cyc = [lon; 360+lon[1]]
DIP_2D_cyc = hcat(DIP_2D, DIP_2D[:,1])

plt = contourf(lon_cyc, lat, DIP_2D_cyc, levels=0:0.2:3.6, transform=ccrs.PlateCarree(), zorder=-1)
colorbar(plt, orientation="horizontal");
title("PO₄ at $(string(round(typeof(1u"m"),grd.depth[iz]))) depth using the OCIM0.1 circulation")
gcf()

using WorldOceanAtlasTools
μDIPobs3D, σ²DIPobs3D = WorldOceanAtlasTools.fit_to_grid(grd, 2018, "phosphate", "annual", "1°", "an")
μDIPobs, σ²DIPobs = μDIPobs3D[iwet], σ²DIPobs3D[iwet]
μx = (μDIPobs, missing)
σ²x = (σ²DIPobs, missing)

ωs = [1.0, 0.0] # the weight for the mismatch (weight of POP = 0)
ωp = 1e-4       # the weight for the parameters prior estimates

v = ustrip.(vector_of_volumes(wet3D, grd))
f   =   generate_objective(ωs, μx, σ²x, v, ωp, mean_obs(p), variance_obs(p))
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

results = optimize(obj, grad, hess, λ, NewtonTrustRegion(), opt)

p

p_optimized = λ2p(results.minimizer)

prob_optimized = SteadyStateProblem(F, ∇ₓF, s, p_optimized)
s_optimized = solve(prob_optimized, CTKAlg()).u

DIPold, _ = state_to_tracers(s, nb, 2)
DIPnew, _ = state_to_tracers(s_optimized, nb, 2)
δDIPold_3D = fill(NaN, size(grd))
δDIPold_3D[iwet] .= 100(DIPold - μDIPobs) ./ μDIPobs
δDIPnew_3D = fill(NaN, size(grd))
δDIPnew_3D[iwet] .= 100(DIPnew - μDIPobs) ./ μDIPobs

δDIPold_2D = δDIPold_3D[:,:,iz]
δDIPnew_2D = δDIPnew_3D[:,:,iz]

δDIPold_cyc = hcat(δDIPold_2D, δDIPold_2D[:,1])
δDIPnew_cyc = hcat(δDIPnew_2D, δDIPnew_2D[:,1])

figure()
δDIPlevels = -20:2:20
ax = subplot(projection = ccrs.EqualEarth(central_longitude=-155.0))
ax.coastlines()
plt2 = contourf(lon_cyc, lat, δDIPold_cyc, cmap="PiYG_r", levels=δDIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="both")
cbar2 = colorbar(plt2, orientation="horizontal", extend="both")
cbar2.set_label("δDIP / DIP [%]")
title("old DIP mismatch at $(string(round(typeof(1u"m"),grd.depth[iz]))) depth using the OCIM0.1 circulation")
gcf()

figure()
ax = subplot(projection = ccrs.EqualEarth(central_longitude=-155.0))
ax.coastlines()
plt3 = contourf(lon_cyc, lat, δDIPnew_cyc, cmap="PiYG_r", levels=δDIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="both")
cbar3 = colorbar(plt3, orientation="horizontal", extend="both")
cbar3.set_label("δDIP / DIP [%]")
title("new DIP mismatch at $(string(round(typeof(1u"m"),grd.depth[iz]))) depth using the OCIM0.1 circulation")
gcf()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

