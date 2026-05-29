#---------------------------------------------------------
# # [Parameter optimisation](@id parameter-optimization)
#---------------------------------------------------------

# In this how-to we fit two parameters of the coupled PO₄–POP model
# (`τ_DIP`, `τ_POP`) to PO₄ observations from
# [WorldOceanAtlasTools.jl](https://github.com/briochemc/WorldOceanAtlasTools.jl)
# on the OCCA circulation. Restricting the fit to two parameters lets
# us also sweep the cost function in a 2D grid around the optimiser
# trajectory and visualise the loss landscape directly.
#
# We use the SciML [Optimization.jl](https://docs.sciml.ai/Optimization/stable/)
# frontend, with [F1Method.jl](https://github.com/briochemc/F1Method.jl)
# supplying the cached steady-state gradient and Hessian, and share an
# UMFPACK factorisation between F1Method's adjoint solve and AIBECS's
# inner Newton step (warm-start trick) so the sweep stays cheap.
#
# The API surface we touch:
#
# - `@flattenable @bounds @logscaled @prior` — declare which parameters
#   are optimised and supply their priors
# - [`p2λ`](@ref) / [`λ2p`](@ref) — transform a parameter struct to/from
#   an unconstrained flat vector
# - [`f_and_∇ₓf`](@ref) — build the volume-weighted mismatch objective
#   and its gradient w.r.t. the state vector x
# - `F1Method` — gradient and Hessian
#   ([Pasquier et al., 2019](https://doi.org/10.5194/gmd-12-3537-2019))
# - `Optimization.jl` — the SciML optimisation frontend, wrapping
#   many optimizers (we will use `NewtonTrustRegion`).
# - [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/) — the SciML
#   linear-solver frontend (we will use the default sparse UMFPACK), and share
#   the factorisation between F1Method and AIBECS.

# ## Re-build the P-model on OCCA

# We reuse the source/sink terms from the P-model tutorial verbatim, just
# swapping the circulation to OCCA so the resolution stays moderate
# (~100k wet boxes × 2 tracers; OCIM2 is ~200k wet boxes for comparison).

using AIBECS
using JLD2 # required by `OCCA.load`
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

# ## Optimisable parameter struct

# We layer the full optimisation metadata stack on top of the tutorial's
# `PmodelParameters`. The leftmost macro consumes the leftmost column
# (canonical example at `src/Parameters.jl`). Only
# `τ_DIP` and `τ_POP` are optimised and marked as `flattenable`.
# So we only pass a 2-element vector `λ` to the optimiser.

import AIBECS: @initial_value, initial_value
import AIBECS: @units, units
import AIBECS: @flattenable, flattenable
import AIBECS: @bounds, bounds
import AIBECS: @logscaled, logscaled
import AIBECS: prior
using Unitful: m, d, s, yr, kyr, Myr, mol, mmol, μmol, μM
using Distributions
using Bijectors    # activates the `AIBECSBijectorsExt` extension (p2λ, λ2p)
using DataFrames   # activates the `AIBECSDataFramesExt` extension (AIBECS.table)

const ∞ = Inf

@initial_value @units @flattenable @bounds @logscaled struct OptParams{U} <: AbstractParameters{U}
    w₀::U      |   5.0  | m/d      | false | (0, ∞) | true
    w′::U      |  0.035 | m/d/m    | false | (0, ∞) | true
    τ_DIP::U   | 100.0  | d        | true  | (0, ∞) | true
    k::U       |  10.0  | μmol/m^3 | false | (0, ∞) | true
    z₀::U      |  80.0  | m        | false | (0, ∞) | false
    τ_POP::U   |  30.0  | d        | true  | (0, ∞) | true
    τ_geo::U   |   1.0  | kyr      | false | (0, ∞) | true
    DIP_geo::U |   2.12 | mmol/m^3 | false | (0, ∞) | false
end

# We apply a penalty to the parameters. For this we
# derive priors from `(lb, ub) = bounds(T, s)` and the initial value, matching
# the convention in `test/parameters.jl` and GNOM: positive half-line gets a
# LogNormal, the real line a wide Normal, and finite-interval parameters get a
# LogitNormal stretched onto `(lb, ub)`.
# Here we will only use the `LogNormal` priors but this gives you an idea
# of how to set up different priors if you want to.

#md # !!! note
#md #     Under the hood AIBECS will perform a change of variables based
#md #     on Turing.jl's Bijectors.jl extension to transform the
#md #     constrained parameter space `(0, ∞)` into an unconstrained one,
#md #     which is what the optimiser operates in. This is why we can use
#md #     `LogNormal` and `LogitNormal` priors even though the optimiser
#md #     only sees unconstrained vectors.

function prior(::Type{T}, s::Symbol) where {T <: OptParams}
    flattenable(T, s) || return nothing
    lb, ub = bounds(T, s)
    if (lb, ub) == (0, ∞)
        return LogNormal(log(initial_value(T, s)), 1.0)
    elseif (lb, ub) == (-∞, ∞)
        return Normal(initial_value(T, s), 10.0)
    else
        m = initial_value(T, s)
        f = (m - lb) / (ub - lb)
        return LocationScale(lb, ub - lb, LogitNormal(log(f / (1 - f)), 1.0))
    end
end

p₀ = OptParams()

# Two fields (`τ_DIP`, `τ_POP`) carry `flattenable = true`; the rest
# are held fixed at their initial value. The metadata table:

AIBECS.table(p₀)

# ## State function in parameter-type form

# Using the `AIBECSFunction(Ts, Gs, nb, ::Type{P})` form lets the SciML
# solvers feed in a flat parameter vector `λ`; AIBECS will call
# `λ2p(OptParams, λ)` internally to reconstruct the struct.

nb = sum(iswet(grd))
F = AIBECSFunction((T_DIP, T_POP), (G_DIP, G_POP), nb, OptParams)
x0 = p₀.DIP_geo * ones(2nb)

# ## Load WOA PO₄ observations

# `WorldOceanAtlasTools.fit_to_grid` returns the mean and variance of
# WOA PO₄ regridded onto our `grd` as 3D fields. There is a subtle
# unit gotcha: the WOA 2018 PO₄ file's `units` attribute is
# `"micromoles_per_kilogram"`, and `WorldOceanAtlasTools.convert_to_SI_unit!`
# does `χ_3D .*= ustrip(upreferred(1.0u"μmol/kg"))`, which gives
# `1.0e-6 mol/kg` — kg is already SI, so `upreferred` does **not**
# convert to `mol/m³`. The returned array is therefore in `mol/kg`
# (mean ≈ 1.6e-6), a factor of ~1000 smaller than the AIBECS model
# state. We close the gap with a nominal seawater density before
# anything downstream consumes these arrays.

using WorldOceanAtlasTools
μDIPobs3D, σ²DIPobs3D = WorldOceanAtlasTools.fit_to_grid(grd, "PO₄")
iwet = iswet(grd)

ρ_factor = ustrip(upreferred(1.035u"kg/L"))   # = 1035.0 kg/m³
μDIPobs  = μDIPobs3D[iwet]  .* ρ_factor       # mol/kg → mol/m³
σ²DIPobs = σ²DIPobs3D[iwet] .* ρ_factor^2     # (mol/kg)² → (mol/m³)²
v = volumevec(grd)

# ## Build the objective `(f, ∇ₓf)`

# We weight the DIP misfit at 1, ignore POP (no direct observations),
# and add a small parameter-prior penalty.

ωs = (1.0, 0.0)
ωp = 1.0e-4
μx  = (μDIPobs,  missing)
σ²x = (σ²DIPobs, missing)
f, ∇ₓf = f_and_∇ₓf(ωs, μx, σ²x, v, ωp, OptParams)

# ## F-1 memory cache + warm-started linear solve

# `F1Method.initialize_mem` solves the steady state once and caches the
# factorised Jacobian; later objective / gradient / Hessian calls reuse it
# whenever the parameter vector has not moved. We route F1Method's linear
# algebra through `LinearSolve.jl` and share its UMFPACK factors with a
# sibling cache that AIBECS's inner Newton step uses
# (`AIBECScache.cacheval = mem.linear_cache.cacheval`) — this makes the
# cost-surface sweep below ~10× faster. `CTKAlg`'s Shamanskii + Armijo
# logic refreshes the factors automatically when they go stale.

using F1Method, Optimization, OptimizationOptimJL, ADTypes
using LinearSolve

solver_kwargs = (; maxItNewton = 15)

λ₀ = p2λ(p₀)
mem = F1Method.initialize_mem(F, ∇ₓf, x0, λ₀, CTKAlg();
    ad = AutoForwardDiff(),
    linsolve = UMFPACKFactorization(),
    solver_kwargs...)

AIBECScache = init(
    LinearProblem(F.jac(x0, λ2p(OptParams, λ₀)), similar(x0)),
    UMFPACKFactorization(),
)
AIBECScache.cacheval = mem.linear_cache.cacheval
AIBECScache.isfresh  = false

solver_kwargs = (; solver_kwargs..., linear_cache = AIBECScache)

optfn = F1Method.optimization_function(f, F, ∇ₓf, mem, CTKAlg(); solver_kwargs...)

# ## Run the optimiser

# `maxiters = 20` is plenty to watch the outer Newton trust-region
# trace drive the gradient norm down by several decades. We attach a
# callback that records the parameter vector at each iteration so the
# diagnostics figure below can overlay the Newton trajectory on the
# cost-surface contourf.

λ_trajectory = Vector{Vector{Float64}}()
push!(λ_trajectory, copy(λ₀))
trajectory_cb = function (state, _)
    push!(λ_trajectory, copy(state.u))
    return false
end

prob = OptimizationProblem(optfn, copy(λ₀))
results = solve(prob, NewtonTrustRegion();
    maxiters = 20, callback = trajectory_cb, show_trace = true)
p_opt = λ2p(OptParams, results.u)

# ## Cost-function sweep around the trajectory

# Evaluate `f(x_ss, λ)` on a 10×10 grid in λ-space centred on the trajectory
# (padded by ±0.2). Each cell solves the steady state from the same
# cold-start `x0`, warm-started by the shared UMFPACK cache from the previous
# section. Takes a couple of minutes.

n_grid = 10
pad    = 0.2
λ1_lo, λ1_hi = extrema(λ[1] for λ in λ_trajectory) .+ (-pad, pad)
λ2_lo, λ2_hi = extrema(λ[2] for λ in λ_trajectory) .+ (-pad, pad)
λ1_vec = range(λ1_lo, λ1_hi; length = n_grid)
λ2_vec = range(λ2_lo, λ2_hi; length = n_grid)

F_grid = Matrix{Float64}(undef, n_grid, n_grid)
for j in 1:n_grid, i in 1:n_grid
    λ = [λ1_vec[i], λ2_vec[j]]
    sol = solve(SteadyStateProblem(F, x0, λ), CTKAlg(); solver_kwargs...)
    F_grid[i, j] = f(sol.u, λ)
end

# ## Diagnostics — setup
#
# We use `CairoMakie` (see also the [Makie how-to](@ref plots-makie)) for plotting.

using CairoMakie
using OceanBasins
using KernelDensity
CairoMakie.activate!(type = "png")

prob_opt = SteadyStateProblem(F, x0, p_opt)
DIPopt, _ = state_to_tracers(solve(prob_opt, CTKAlg(); solver_kwargs...).u, grd)

DIPopt_μM = ustrip.(μM, DIPopt  .* (mol/m^3))
μWOA_μM   = ustrip.(μM, μDIPobs .* (mol/m^3))
diff_μM   = DIPopt_μM .- μWOA_μM

depth = ustrip.(grd.depth)
lat   = ustrip.(grd.lat)
lon   = ustrip.(grd.lon)
maxdepth_round = ceil(maximum(depth) / 1000) * 1000

dip_levels  = 0:0.25:3.0
dip_cmap    = cgrad(:viridis, length(dip_levels) - 1; categorical = true)
diff_levels = -1.0:0.2:1.0
diff_cmap   = cgrad(:RdBu, length(diff_levels) - 1; categorical = true, rev = true)

latticks = (-90:30:90, ["90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"])
lonticks = (0:60:360, ["0°", "60°E", "120°E", "180°", "120°W", "60°W", "0°"])
gcdepthticks = (0:1000:maxdepth_round, string.(Int.(0:1000:maxdepth_round)))
plus_minus_ticks = (-1:0.5:1, ["−1", "−0.5", "0", "+0.5", "+1"])

# Basin masks, each extended to the pole via its Southern Ocean sector.
OCEANS = oceanpolygons()
lat_grd, lon_grd = latvec(grd), lonvec(grd)
mPAC = ispacific2(lat_grd, lon_grd, OCEANS)  .| isSOpacific(lat_grd, lon_grd, OCEANS)
mATL = isatlantic2(lat_grd, lon_grd, OCEANS) .| isSOatlantic(lat_grd, lon_grd, OCEANS)
mIND = isindian2(lat_grd, lon_grd, OCEANS)   .| isSOindian(lat_grd, lon_grd, OCEANS)
basin_specs = [
    (name = "Pacific",  mask = mPAC, color = :dodgerblue),
    (name = "Atlantic", mask = mATL, color = :firebrick),
    (name = "Indian",   mask = mIND, color = :darkorange),
]

# ## Optimised vs WOA scatter
#
# Each wet box at `(WOA, optimised)`, coloured by its mass-weighted density
# quantile under a volume-weighted bivariate KDE. Walking the sorted
# cumulative density once gives both the iso-mass contour levels and the
# per-point colours.

bandwidth = 0.25 .* KernelDensity.default_bandwidth((μWOA_μM, DIPopt_μM))
d_kde = kde((μWOA_μM, DIPopt_μM); bandwidth, weights = v ./ sum(v))
ik = InterpKDE(d_kde)

dx, dy = step(d_kde.x), step(d_kde.y)
ρ_asc = sort(vec(d_kde.density))
cum_asc = cumsum(ρ_asc) .* (dx * dy)
cum_asc ./= cum_asc[end]
mass_qs = 0.1:0.1:0.9
levels = [ρ_asc[searchsortedfirst(cum_asc, q)] for q in mass_qs]
point_q = [(i = searchsortedlast(ρ_asc, d);
            i == 0 ? 0.0 : cum_asc[i]) for d in ik.itp.(μWOA_μM, DIPopt_μM)]

fig = Figure()
ax = Axis(fig[1, 1];
    xlabel = "WOA PO₄ (μM)", ylabel = "optimised PO₄ (μM)",
    xgridvisible = false, ygridvisible = false,
    aspect = DataAspect(), limits = (0, 3, 0, 3))
sc = scatter!(ax, μWOA_μM, DIPopt_μM;
    markersize = 4, color = point_q, colormap = :plasma, colorrange = (0, 1))
contour!(ax, d_kde.x, d_kde.y, d_kde.density;
    levels = levels, color = (:black, 0.5), linewidth = 1.2)
lines!(ax, [0, 3], [0, 3]; color = :red, linestyle = :dash)
Colorbar(fig[1, 2], sc; label = "density quantile", ticks = mass_qs)
fig

# ## Per-basin depth profiles — model vs WOA

fig = Figure(size = (900, 350))
for (c, basin) in enumerate(basin_specs)
    ax = Axis(fig[1, c];
        xgridvisible = false, ygridvisible = false,
        xlabel = "PO₄ (μM)",
        ylabel = c == 1 ? "depth (m)" : "",
        yreversed = true, yticks = gcdepthticks,
        title = basin.name,
        limits = (nothing, (0, maxdepth_round)))
    prof_model = AIBECS.horizontalmean(DIPopt_μM, grd, basin.mask)
    prof_obs   = AIBECS.horizontalmean(μWOA_μM,   grd, basin.mask)
    lines!(ax, prof_model, depth;
        color = basin.color, linewidth = 2, label = "model")
    lines!(ax, prof_obs, depth;
        color = basin.color, linewidth = 2, linestyle = :dash, label = "WOA")
    c == 1 && axislegend(ax; position = :rb, framevisible = false)
    c > 1 && hideydecorations!(ax;
        ticklabels = true, label = true, ticks = false, grid = false)
end
fig

# ## Per-basin zonal means — model / WOA / model − WOA

row_specs = [
    (label = "model",         field = DIPopt_μM, levels = dip_levels,  cmap = dip_cmap),
    (label = "WOA",           field = μWOA_μM,   levels = dip_levels,  cmap = dip_cmap),
    (label = "model − WOA",   field = diff_μM,   levels = diff_levels, cmap = diff_cmap),
]
zm_grid = [
    [AIBECS.zonalmean(spec.field, grd, basin.mask) for basin in basin_specs]
    for spec in row_specs
]
function wetlatrange(zm, pad = 5)
    haswet = [any(!isnan, view(zm, i, :)) for i in eachindex(lat)]
    wetidx = findall(haswet)
    isempty(wetidx) && return (-90.0, 90.0)
    return (max(-90, lat[first(wetidx)] - pad),
            min( 90, lat[last(wetidx)]  + pad))
end
xlim_per_basin = [wetlatrange(zm_grid[1][c]) for c in eachindex(basin_specs)]

## Documenter `@example` blocks run at module scope, where bare
## `co_v = co` inside a `for` would shadow the outer binding; `Ref` sidesteps it.
co_v = Ref{Any}()
co_d = Ref{Any}()

fig = Figure(size = (1400, 900))
for (r, spec) in enumerate(row_specs), (c, basin) in enumerate(basin_specs)
    xlim = xlim_per_basin[c]
    ax = Axis(fig[r, c];
        backgroundcolor = :lightgray,
        xgridvisible = false, ygridvisible = false,
        xticksmirrored = true, yticksmirrored = true,
        xticks = latticks, yticks = gcdepthticks,
        yreversed = true,
        ylabel = c == 1 ? "depth (m)" : "",
        title = r == 1 ? basin.name : "",
        limits = (xlim[1], xlim[2], 0, maxdepth_round))
    co = contourf!(ax, lat, depth, zm_grid[r][c];
        levels = spec.levels, colormap = spec.cmap, nan_color = :lightgray,
        extendlow = spec.cmap[1], extendhigh = spec.cmap[end])
    r == 1 && c == 1 && (co_v[] = co)
    r == 3 && c == 1 && (co_d[] = co)
    c > 1 && hideydecorations!(ax; ticklabels = true, label = true, ticks = false, grid = false)
    r < 3 && hidexdecorations!(ax; ticklabels = true, label = true, ticks = false, grid = false)
end
Label(fig[1, 0], "model";        rotation = π/2, font = :bold, fontsize = 16, tellheight = false)
Label(fig[2, 0], "WOA";          rotation = π/2, font = :bold, fontsize = 16, tellheight = false)
Label(fig[3, 0], "model − WOA";  rotation = π/2, font = :bold, fontsize = 16, tellheight = false)
Colorbar(fig[1:2, length(basin_specs) + 1], co_v[];
    label = "PO₄ (μM)", ticks = 0:1:3)
Colorbar(fig[3, length(basin_specs) + 1], co_d[];
    label = "model − WOA (μM)", ticks = plus_minus_ticks)
colsize!(fig.layout, 0, Auto(0.15))
for (c, xlim) in enumerate(xlim_per_basin)
    colsize!(fig.layout, c, Auto(xlim[2] - xlim[1]))
end
colsize!(fig.layout, length(basin_specs) + 1, Auto(0.25))
fig

# ## Horizontal slices at 1000 m — model / WOA / diff

sl_opt  = AIBECS.horizontalslice(DIPopt_μM, grd; depth = 1000)
sl_obs  = AIBECS.horizontalslice(μWOA_μM,   grd; depth = 1000)
sl_diff = sl_opt .- sl_obs

function map_axis(grid_pos; title = "")
    Axis(grid_pos;
        backgroundcolor = :lightgray,
        xgridvisible = false, ygridvisible = false,
        xticks = lonticks, yticks = latticks,
        title = title, limits = (0, 360, -90, 90))
end

fig = Figure(size = (1400, 500))
ax_m = map_axis(fig[1, 1]; title = "optimised")
co_m = contourf!(ax_m, lon, lat, sl_opt';
    levels = dip_levels, colormap = dip_cmap, nan_color = :lightgray,
    extendlow = dip_cmap[1], extendhigh = dip_cmap[end])
ax_o = map_axis(fig[1, 2]; title = "WOA")
contourf!(ax_o, lon, lat, sl_obs';
    levels = dip_levels, colormap = dip_cmap, nan_color = :lightgray,
    extendlow = dip_cmap[1], extendhigh = dip_cmap[end])
ax_d = map_axis(fig[1, 3]; title = "optimised − WOA")
co_dd = contourf!(ax_d, lon, lat, sl_diff';
    levels = diff_levels, colormap = diff_cmap, nan_color = :lightgray,
    extendlow = diff_cmap[1], extendhigh = diff_cmap[end])
linkyaxes!(ax_m, ax_o, ax_d)
hideydecorations!(ax_o; ticklabels = true, label = false, ticks = false, grid = false)
hideydecorations!(ax_d; ticklabels = true, label = false, ticks = false, grid = false)
cb_v = Colorbar(fig[2, 1:2], co_m; vertical = false, flipaxis = false,
    label = "PO₄ @ 1000 m (μM)", ticks = 0:1:3)
cb_v.width = Relative(0.7)
cb_diff = Colorbar(fig[2, 3], co_dd; vertical = false, flipaxis = false,
    label = "optimised − WOA (μM)", ticks = plus_minus_ticks)
cb_diff.width = Relative(0.85)
fig

# ## Side-by-side parameter tables, initial vs optimised

AIBECS.table(p₀)

#-

AIBECS.table(p_opt)

# ## Extra: cost-surface sweep + Newton trajectory
#
# Mostly diagnostic. Once you've converged, the optimal parameters are the
# point — but the sweep + trajectory is useful when something goes wrong.
# Convert the λ-space sweep grid and trajectory to physical (τ_DIP, τ_POP)
# via `λ2p`, then a single `contourf!` + `scatterlines!` does it.

τDIP_grid = [λ2p(OptParams, [λ1, λ₀[2]]).τ_DIP for λ1 in λ1_vec]
τPOP_grid = [λ2p(OptParams, [λ₀[1], λ2]).τ_POP for λ2 in λ2_vec]
τDIP_traj = [λ2p(OptParams, λ).τ_DIP for λ in λ_trajectory]
τPOP_traj = [λ2p(OptParams, λ).τ_POP for λ in λ_trajectory]

fig = Figure()
ax = Axis(fig[1, 1];
    xlabel = "τ_DIP (d)", ylabel = "τ_POP (d)",
    xscale = log10, yscale = log10,
    title = "√f(p) cost surface + Newton trajectory")
co = contourf!(ax, τDIP_grid, τPOP_grid, sqrt.(F_grid);
    levels = 20, colormap = cgrad(:BrBG; rev = true))
scatterlines!(ax, τDIP_traj, τPOP_traj;
    color = :red, linewidth = 2, linestyle = :dash, marker = :circle, markersize = 4)
scatter!(ax, [τDIP_traj[end]], [τPOP_traj[end]];
    marker = :star5, markersize = 22, color = :yellow,
    strokecolor = :black, strokewidth = 1.5)
annotation!(ax, 70, -40, τDIP_traj[1], τPOP_traj[1];
    text = "first iterate", path = Ann.Paths.Corner())
annotation!(ax, -70, +40, τDIP_traj[end], τPOP_traj[end];
    text = "optimal solution", path = Ann.Paths.Corner(), color = :white)
Colorbar(fig[1, 2], co; label = "√f(p)")
fig

# ## Where to next

# This example only fits two parameters of a single-tracer mismatch
# objective on a moderate-resolution circulation; the multi-tracer,
# multi-observation fit lives in [GNOM](https://github.com/MTEL-USC/GNOM).
# See the [roadmap](@ref roadmap) for alternative parameter-estimation
# backends (Turing, ImplicitDifferentiation, SciMLSensitivity), the
# table-form `f_and_∇ₓf(ωs, ωp, grd, obs, ::Type{P})` for fitting against
# `WorldOceanAtlasTools.observations(...)` tables, and other directions
# where contributions are welcome.
