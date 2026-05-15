# Cost-function sweep + optimizer-trajectory overlay
# ====================================================
#
# Exploratory companion to `docs/lit/howtos/4_parameter_optimization.jl`.
# Not part of the rendered docs — run manually:
#
#     julia --project=docs docs/experiments/cost_sweep.jl
#
# What it does:
#  1. Same P-model on OCCA as the how-to, but with only two flattenable
#     parameters (`τ_DIP`, `τ_POP`) so the cost is a function of a 2D λ.
#  2. Sweeps a 10×10 grid of λ values centred on `p2λ(p₀)` (one fresh
#     nonlinear solve per grid point — ~100 solves, ~10–15 min on OCCA).
#  3. Runs `NewtonTrustRegion` from `p₀`, recording the λ at every
#     optimiser step.
#  4. Saves a `contourf` of `f(λ)` with the optimiser trajectory and
#     start/end points overlaid to `docs/experiments/cost_sweep.png`.

using AIBECS
using JLD2
using Unitful: m, d, s, yr, kyr, Myr, mol, mmol, μmol, μM
using Distributions
using Bijectors
using DataFrames
using WorldOceanAtlasTools
using F1Method
using Optimization
using OptimizationOptimJL
using ADTypes
using Plots
using LinearAlgebra: norm
using LinearSolve
using OceanBasins

# ── Model assembly (mirrors the how-to) ──────────────────────────────────────
grd, T_OCIM = OCCA.load()
T_DIP(p) = T_OCIM
T_POP(p) = transportoperator(grd, z -> w(z, p))

function w(z, p)
    @unpack w₀, w′ = p
    return @. w₀ + w′ * z
end

z = depthvec(grd)
function U(x, p)
    @unpack τ_DIP, k, z₀ = p
    return @. x / τ_DIP * x / (x + k) * (z ≤ z₀) * (x ≥ 0)
end

function R(x, p)
    @unpack τ_POP = p
    return x / τ_POP
end

function G_DIP(DIP, POP, p)
    @unpack DIP_geo, τ_geo = p
    return @. -$U(DIP, p) + $R(POP, p) + (DIP_geo - DIP) / τ_geo
end
function G_POP(DIP, POP, p)
    @unpack τ_geo = p
    return @. $U(DIP, p) - $R(POP, p) - POP / τ_geo
end

import AIBECS: @initial_value, initial_value
import AIBECS: @units, units
import AIBECS: @flattenable, flattenable
import AIBECS: @bounds, bounds
import AIBECS: @logscaled, logscaled
import AIBECS: prior

const ∞ = Inf

# Two-parameter variant: only `τ_DIP` and `τ_POP` carry flattenable=true.
# Both are positive half-line bounded so `p2λ` is just elementwise `log`,
# which makes the sweep axes physically interpretable.
@initial_value @units @flattenable @bounds @logscaled struct OptParams2{U} <: AbstractParameters{U}
    w₀::U      |   5.0  | m/d      | false | (0, ∞) | true
    w′::U      |  0.035 | m/d/m    | false | (0, ∞) | true
    τ_DIP::U   |  100.0 | d        | true  | (0, ∞) | true
    k::U       |  10.0  | μmol/m^3 | false | (0, ∞) | true
    z₀::U      |  80.0  | m        | false | (0, ∞) | false
    τ_POP::U   |  30.0  | d        | true  | (0, ∞) | true
    τ_geo::U   |   1.0  | kyr      | false | (0, ∞) | true
    DIP_geo::U |   2.12 | mmol/m^3 | false | (0, ∞) | false
end

function prior(::Type{T}, s::Symbol) where {T <: OptParams2}
    flattenable(T, s) || return nothing
    lb, ub = bounds(T, s)
    if (lb, ub) == (0, ∞)
        return LogNormal(log(initial_value(T, s)), 1.0)
    elseif (lb, ub) == (-∞, ∞)
        return Normal(initial_value(T, s), 10.0)
    else
        m = initial_value(T, s)
        f_ = (m - lb) / (ub - lb)
        return LocationScale(lb, ub - lb, LogitNormal(log(f_ / (1 - f_)), 1.0))
    end
end

p₀ = OptParams2()
nb = sum(iswet(grd))
F = AIBECSFunction((T_DIP, T_POP), (G_DIP, G_POP), nb, OptParams2)
@unpack DIP_geo = p₀
x0 = DIP_geo * ones(2nb)

# Sanity-check mass conservation for both transport operators
# (`norm(v) / norm(Tᵀv) > 1 Myr` is the test/setup.jl criterion).
using LinearAlgebra: norm
let v_vec = volumevec(grd), Myr_in_s = ustrip(upreferred(1u"Myr"))
    T_DIP_mat = T_DIP(p₀)
    T_POP_mat = T_POP(p₀)
    @assert norm(v_vec) / norm(T_DIP_mat' * v_vec) > Myr_in_s "T_DIP does not conserve mass"
    @assert norm(v_vec) / norm(T_POP_mat' * v_vec) > Myr_in_s "T_POP does not conserve mass"
    @info "Mass-conservation timescales (should exceed 1 Myr)" T_DIP_Myr = norm(v_vec) / norm(T_DIP_mat' * v_vec) / Myr_in_s T_POP_Myr = norm(v_vec) / norm(T_POP_mat' * v_vec) / Myr_in_s
end

μDIPobs3D, σ²DIPobs3D = WorldOceanAtlasTools.fit_to_grid(grd, "PO₄")
iwet = iswet(grd)

# WOA 2018 PO₄ file's `units` attribute is "micromoles_per_kilogram".
# `WorldOceanAtlasTools.convert_to_SI_unit!` does
# `χ_3D .*= ustrip(upreferred(1.0u"μmol/kg"))`, which gives `1.0e-6 mol/kg`
# — kg is already SI, so `upreferred` does NOT convert to mol/m³. The
# returned `μDIPobs` is therefore in mol/kg (mean ≈ 1.6e-6), a factor of
# ~1000 smaller than the AIBECS model state in mol/m³. We close the gap
# with a nominal seawater density before everything downstream
# (cost function, plotting) consumes these arrays.
ρ_factor = ustrip(upreferred(1.035u"kg/L"))   # = 1035.0 kg/m^3
μDIPobs  = μDIPobs3D[iwet]  .* ρ_factor       # mol/kg → mol/m³
σ²DIPobs = σ²DIPobs3D[iwet] .* ρ_factor^2     # (mol/kg)² → (mol/m³)²
v = volumevec(grd)

ωs = (1.0, 0.0)
ωp = 1.0e-4
μx  = (μDIPobs,  missing)
σ²x = (σ²DIPobs, missing)
f, ∇ₓf = f_and_∇ₓf(ωs, μx, σ²x, v, ωp, OptParams2)

λ₀ = p2λ(p₀)
@assert length(λ₀) == 2 "Expected 2 flattenable parameters; got $(length(λ₀))"
@info "Initial λ" λ₀

solver_kwargs = (; maxItNewton = 15)

# Route F1Method's internal linear algebra through LinearSolve.jl so the
# factorised Jacobian lives in a `LinearSolve.LinearCache` we can share
# with AIBECS's inner Newton solver below.
mem = F1Method.initialize_mem(F, ∇ₓf, x0, λ₀, CTKAlg();
                              ad = AutoForwardDiff(),
                              linsolve = UMFPACKFactorization(),
                              solver_kwargs...)

# Warm-start cache for AIBECS's Newton step. F1Method's `mem.linear_cache`
# carries a matrix RHS (the parameter Jacobian `-∇ₚF`); AIBECS needs a
# vector RHS, so we keep a sibling cache and share the UMFPACK factors via
# `cacheval`. Built once at λ₀ and threaded through *every* downstream
# CTKAlg solve (sweep + optimisation). When the warm-start gets too stale
# CTKAlg's Shamanskii ratio check + Armijo line search refresh the factors
# automatically, costing at most one extra Newton iteration.
AIBECScache = init(
    LinearProblem(F.jac(x0, λ2p(OptParams2, λ₀)), similar(x0)),
    UMFPACKFactorization(),
)
AIBECScache.cacheval = mem.linear_cache.cacheval
AIBECScache.isfresh  = false

solver_kwargs = (; solver_kwargs..., linear_cache = AIBECScache)

# ── Optimization with trajectory recording ──────────────────────────────────
# Run the optimizer first so we can size the sweep grid around the trajectory,
# guaranteeing the optimizer path is rendered inside the `contourf`.
λ_trajectory = Vector{Vector{Float64}}()
push!(λ_trajectory, copy(λ₀))

trajectory_cb = function (opt_state, _)
    push!(λ_trajectory, copy(opt_state.u))
    return false
end

optfn = F1Method.optimization_function(f, F, ∇ₓf, mem, CTKAlg(); solver_kwargs...)
prob  = OptimizationProblem(optfn, copy(λ₀))
results = solve(prob, NewtonTrustRegion();
                maxiters = 20,
                show_trace = false,
                callback = trajectory_cb)
p_opt = λ2p(OptParams2, results.u)

@info "Optimizer finished" retcode = results.retcode n_steps = length(λ_trajectory) p_opt

# ── Sweep in λ-space, sized around the optimizer trajectory ─────────────────
# Bounds come from `extrema(λ_trajectory)` padded by ±0.2 in log-space, so
# the trajectory always sits inside the contour. Every cell is solved
# independently from the same cold-start `x0`. The retcode is ignored —
# even when Newton hits `maxItNewton` without fully converging, the residual
# is small enough that `f(u, λ)` is a fine approximation of the true cost
# for a coarse contour.

n_grid = 10
pad    = 0.2
λ1_lo, λ1_hi = extrema(λ[1] for λ in λ_trajectory) .+ (-pad, pad)
λ2_lo, λ2_hi = extrema(λ[2] for λ in λ_trajectory) .+ (-pad, pad)
λ1_vec = range(λ1_lo, λ1_hi; length = n_grid)
λ2_vec = range(λ2_lo, λ2_hi; length = n_grid)

F_grid = Matrix{Float64}(undef, n_grid, n_grid)
@info "Starting $(n_grid)×$(n_grid) sweep ($(n_grid * n_grid) solves)..."
t_sweep = @elapsed for j in 1:n_grid, i in 1:n_grid
    λ = [λ1_vec[i], λ2_vec[j]]
    local p_struct = λ2p(OptParams2, λ)
    @info " τ_DIP = $(round(p_struct.τ_DIP; digits=2)), τ_POP = $(round(p_struct.τ_POP; digits=2))"

    sol = solve(SteadyStateProblem(F, x0, λ), CTKAlg();
                preprint = "    ", solver_kwargs...)
    F_grid[i, j] = f(sol.u, λ)
end
@info "Sweep complete" t_sweep_s = round(t_sweep; digits = 1)

# ── Plot: contourf of f in p-space (log10 axes) + optimizer trajectory ─────
# F_grid is indexed as F_grid[i, j] for (λ1_vec[i], λ2_vec[j]).
# Plots.contourf expects z[row, col] = z(x[col], y[row]), so transpose.
# Use `λ2p` to convert λ-space coordinates to the parameter struct, then
# pull out τ_DIP and τ_POP — this stays correct regardless of which
# flattenable transform (log / logit / identity) each parameter uses.

τDIP_vec  = [λ2p(OptParams2, [λ1, λ₀[2]]).τ_DIP for λ1 in λ1_vec]
τPOP_vec  = [λ2p(OptParams2, [λ₀[1], λ2]).τ_POP for λ2 in λ2_vec]
τDIP_traj = [λ2p(OptParams2, λ).τ_DIP for λ in λ_trajectory]
τPOP_traj = [λ2p(OptParams2, λ).τ_POP for λ in λ_trajectory]

plt = contourf(τDIP_vec, τPOP_vec, sqrt.(F_grid');
               xlabel = "τ_DIP [d]",
               ylabel = "τ_POP [d]",
               title  = "√f(p) sweep + NewtonTrustRegion trajectory",
               levels = 20,
               c = :viridis,
               colorbar_title = "√f")
plot!(xscale = :log10, yscale = :log10, minorgrid = true)
plot!(plt, τDIP_traj, τPOP_traj;
      marker     = :circle,
      markersize = 4,
      linewidth  = 2,
      color      = :red,
      label      = "trajectory")
scatter!(plt, [τDIP_traj[1]], [τPOP_traj[1]];
         marker = :star5, markersize = 10, color = :white,
         markerstrokecolor = :black, label = "p₀")
scatter!(plt, [λ2p(OptParams2, results.u).τ_DIP],
              [λ2p(OptParams2, results.u).τ_POP];
         marker = :star5, markersize = 10, color = :red,
         markerstrokecolor = :black, label = "p_opt")

out_png = joinpath(@__DIR__, "cost_sweep.png")
savefig(plt, out_png)
@info "Saved figure" out_png

# ── DIP diagnostics at p_opt ─────────────────────────────────────────────────
# Solve once more at the optimum (warm-started via AIBECScache, same as the
# inner solves), then build three figures comparing the optimised DIP field
# to the WOA observations: horizontal slice, basin-mean profiles, and zonal
# averages per basin.

prob_opt = SteadyStateProblem(F, x0, p_opt)
DIP_opt, _ = state_to_tracers(solve(prob_opt, CTKAlg(); solver_kwargs...).u, grd)

DIPmodel = DIP_opt * (mol / m^3) .|> mmol/m^3
DIPobs   = μDIPobs * (mol / m^3) .|> mmol/m^3
DIPdiff  = DIPmodel - DIPobs

OCEANS = oceanpolygons()
lat_grd, lon_grd = latvec(grd), lonvec(grd)

# Profile basins: Antarctic excluded from the three main basins and plotted
# separately, so each profile is a clean basin-mean without Southern Ocean
# water masses muddying the deep signal.
profile_basins = [
    "Pacific"   => ispacific2(lat_grd, lon_grd, OCEANS),
    "Atlantic"  => isatlantic2(lat_grd, lon_grd, OCEANS),
    "Indian"    => isindian2(lat_grd, lon_grd, OCEANS),
    "Antarctic" => isantarctic(lat_grd, lon_grd, OCEANS),
]

# Zonal-average basins: each basin extended to the pole via its Southern
# Ocean sector, so the zonal sections show every latitude band continuously.
zonal_basins = [
    "Pacific"  => ispacific2(lat_grd, lon_grd, OCEANS)  .| isSOpacific(lat_grd, lon_grd, OCEANS),
    "Atlantic" => isatlantic2(lat_grd, lon_grd, OCEANS) .| isSOatlantic(lat_grd, lon_grd, OCEANS),
    "Indian"   => isindian2(lat_grd, lon_grd, OCEANS)   .| isSOindian(lat_grd, lon_grd, OCEANS),
]

# ── Horizontal slices @ 1000 m: model | obs | diff ──────────────────────────
# Pass numeric clims (Plots strips units when comparing to Unitful tracers).
# The diff clim is hand-set to a tight range (mmol/m³) so moderate-magnitude
# errors aren't washed out by the few extreme cells; the data-driven
# `(-diff_max, diff_max)` saturates the bulk of the field.
slice_depth = 1000m
field_max   = ustrip(maximum(DIPobs))
field_clim  = (0.0, field_max)
diff_clim   = (-1.0, 1.0)
slice_plt = plot(
    plothorizontalslice(DIPmodel, grd; depth = slice_depth, color = :viridis,
                        clim = field_clim, title = "DIP model @ 1000 m"),
    plothorizontalslice(DIPobs,   grd; depth = slice_depth, color = :viridis,
                        clim = field_clim, title = "DIP WOA @ 1000 m"),
    plothorizontalslice(DIPdiff,  grd; depth = slice_depth, color = :balance,
                        clim = diff_clim,  title = "model − WOA"),
    ; layout = (3, 1), size = (900, 1200)
)
slice_png = joinpath(@__DIR__, "cost_sweep_DIP_slice.png")
savefig(slice_plt, slice_png)
@info "Saved DIP slice figure" slice_png

# ── Basin-mean profiles: model vs WOA ───────────────────────────────────────
profile_panels = [
    let p_ = plothorizontalmean(DIPmodel, grd; mask = m, label = "model",
                                color = :black, lw = 2, title = name,
                                xlims = (0, Inf))
        plothorizontalmean!(p_, DIPobs, grd; mask = m, label = "WOA",
                            color = :red, linestyle = :dash, lw = 2)
        p_
    end
    for (name, m) in profile_basins
]
profile_plt = plot(profile_panels...; layout = (2, 2), size = (1000, 800))
profile_png = joinpath(@__DIR__, "cost_sweep_DIP_profiles.png")
savefig(profile_plt, profile_png)
@info "Saved DIP basin-profile figure" profile_png

# ── Zonal averages per basin: model | obs | diff ────────────────────────────
zonal_panels = []
for (name, m) in zonal_basins
    push!(zonal_panels,
        plotzonalaverage(DIPmodel, grd; mask = m, color = :viridis,
                         clim = field_clim, title = "$name model"),
        plotzonalaverage(DIPobs,   grd; mask = m, color = :viridis,
                         clim = field_clim, title = "$name WOA"),
        plotzonalaverage(DIPdiff,  grd; mask = m, color = :balance,
                         clim = diff_clim, title = "$name model − WOA"),
    )
end
zonal_plt = plot(zonal_panels...; layout = (length(zonal_basins), 3), size = (1500, 1600))
zonal_png = joinpath(@__DIR__, "cost_sweep_DIP_zonal.png")
savefig(zonal_plt, zonal_png)
@info "Saved DIP zonal-average figure" zonal_png
