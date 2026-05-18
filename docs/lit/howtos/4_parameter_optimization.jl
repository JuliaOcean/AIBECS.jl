# #---------------------------------------------------------
# # # [Parameter optimisation](@id parameter-optimization)
# #---------------------------------------------------------

# # In this how-to we fit two parameters of the coupled PO₄–POP model
# # (`τ_DIP`, `τ_POP`) to PO₄ observations from
# # [WorldOceanAtlasTools.jl](https://github.com/briochemc/WorldOceanAtlasTools.jl)
# # on the OCCA circulation. Restricting the fit to two parameters lets
# # us also sweep the cost function in a 2D grid around the optimiser
# # trajectory and visualise the loss landscape directly.
# #
# # We use the SciML [Optimization.jl][1] frontend, with [F1Method.jl][2]
# # supplying the cached steady-state gradient and Hessian, and share an
# # UMFPACK factorisation between F1Method's adjoint solve and AIBECS's
# # inner Newton step (warm-start trick) so the sweep stays cheap.
# #
# # The API surface we touch:
# #
# # - `@flattenable @bounds @logscaled @prior` — declare which parameters
# #   are optimised and supply their priors;
# # - [`p2λ`](@ref) / [`λ2p`](@ref) — transform a parameter struct to/from
# #   an unconstrained flat vector;
# # - [`f_and_∇ₓf`](@ref) — build the volume-weighted mismatch objective
# #   and its analytical state Jacobian;
# # - `F1Method` — cached steady-state solve plus adjoint gradient and
# #   Hessian ([Pasquier et al., 2019](https://doi.org/10.5194/gmd-12-3537-2019));
# # - `Optimization.jl` — the SciML optimisation frontend, in the same way
# #   `LinearSolve.jl` is the SciML linear-solver frontend used elsewhere
# #   in AIBECS;
# # - [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/) — we
# #   share an UMFPACK factorisation between F1Method's adjoint solve and
# #   AIBECS's inner Newton step.
# #
# # See the [coupled PO₄–POP tutorial](@ref P-model) for the model
# # equations and [GNOM](https://github.com/MTEL-USC/GNOM) for the
# # converged multi-tracer / isotope fit this example is a thin slice of.
# #
# # [1]: https://docs.sciml.ai/Optimization/stable/
# # [2]: https://github.com/briochemc/F1Method.jl

# # ## Re-build the P-model on OCCA

# # We reuse the source/sink terms from the P-model tutorial verbatim, just
# # swapping the circulation to OCCA so the resolution stays moderate
# # (~170k wet boxes × 2 tracers).

# using AIBECS
# using JLD2 # required by `OCCA.load`
# grd, T_OCIM = OCCA.load()
# T_DIP(p) = T_OCIM
# T_POP(p) = transportoperator(grd, z -> w(z, p))

# function w(z, p)
#     @unpack w₀, w′ = p
#     return @. w₀ + w′ * z
# end

# z = depthvec(grd)
# function U(x, p)
#     @unpack τ_DIP, k, z₀ = p
#     return @. x / τ_DIP * x / (x + k) * (z ≤ z₀) * (x ≥ 0)
# end

# function R(x, p)
#     @unpack τ_POP = p
#     return x / τ_POP
# end

# function G_DIP(DIP, POP, p)
#     @unpack DIP_geo, τ_geo = p
#     return @. -$U(DIP, p) + $R(POP, p) + (DIP_geo - DIP) / τ_geo
# end
# function G_POP(DIP, POP, p)
#     @unpack τ_geo = p
#     return @. $U(DIP, p) - $R(POP, p) - POP / τ_geo
# end

# # ## Optimisable parameter struct

# # We layer the full optimisation metadata stack on top of the tutorial's
# # `PmodelParameters`. The leftmost macro consumes the leftmost column
# # (canonical example at `src/Parameters.jl`); the columns therefore are
# # `initial_value | units | flattenable | bounds | logscaled`. Only two
# # fields — `τ_DIP` and `τ_POP` — carry `flattenable = true`, so the
# # unconstrained vector `λ` has length 2 and the cost-surface sweep
# # below is a 2D contour plot.

# import AIBECS: @initial_value, initial_value
# import AIBECS: @units, units
# import AIBECS: @flattenable, flattenable
# import AIBECS: @bounds, bounds
# import AIBECS: @logscaled, logscaled
# import AIBECS: prior
# using Unitful: m, d, s, yr, kyr, Myr, mol, mmol, μmol, μM
# using Distributions
# using Bijectors    # activates the `AIBECSBijectorsExt` extension (p2λ, λ2p)
# using DataFrames   # activates the `AIBECSDataFramesExt` extension (AIBECS.table)

# const ∞ = Inf

# @initial_value @units @flattenable @bounds @logscaled struct OptParams{U} <: AbstractParameters{U}
#     w₀::U      |   5.0  | m/d      | false | (0, ∞) | true
#     w′::U      |  0.035 | m/d/m    | false | (0, ∞) | true
#     τ_DIP::U   | 100.0  | d        | true  | (0, ∞) | true
#     k::U       |  10.0  | μmol/m^3 | false | (0, ∞) | true
#     z₀::U      |  80.0  | m        | false | (0, ∞) | false
#     τ_POP::U   |  30.0  | d        | true  | (0, ∞) | true
#     τ_geo::U   |   1.0  | kyr      | false | (0, ∞) | true
#     DIP_geo::U |   2.12 | mmol/m^3 | false | (0, ∞) | false
# end

# # Derive priors from `(lb, ub) = bounds(T, s)` and the initial value, matching
# # the convention in `test/parameters.jl` and GNOM: positive half-line gets a
# # LogNormal, the real line a wide Normal, and finite-interval parameters get a
# # LogitNormal stretched onto `(lb, ub)`.

# function prior(::Type{T}, s::Symbol) where {T <: OptParams}
#     flattenable(T, s) || return nothing
#     lb, ub = bounds(T, s)
#     if (lb, ub) == (0, ∞)
#         return LogNormal(log(initial_value(T, s)), 1.0)
#     elseif (lb, ub) == (-∞, ∞)
#         return Normal(initial_value(T, s), 10.0)
#     else
#         m = initial_value(T, s)
#         f = (m - lb) / (ub - lb)
#         return LocationScale(lb, ub - lb, LogitNormal(log(f / (1 - f)), 1.0))
#     end
# end

# p₀ = OptParams()

# # Two fields (`τ_DIP`, `τ_POP`) carry `flattenable = true`; the rest
# # are held fixed at their initial value. The metadata table:

# AIBECS.table(p₀)

# # ## State function in parameter-type form

# # Using the `AIBECSFunction(Ts, Gs, nb, ::Type{P})` form lets the SciML
# # solvers feed in a flat unconstrained vector `λ`; AIBECS will call
# # `λ2p(OptParams, λ)` internally to reconstruct the struct.

# nb = sum(iswet(grd))
# F = AIBECSFunction((T_DIP, T_POP), (G_DIP, G_POP), nb, OptParams)
# x0 = p₀.DIP_geo * ones(2nb)

# # ## Load WOA PO₄ observations

# # `WorldOceanAtlasTools.fit_to_grid` returns the mean and variance of
# # WOA PO₄ regridded onto our `grd` as 3D fields. There is a subtle
# # unit gotcha: the WOA 2018 PO₄ file's `units` attribute is
# # `"micromoles_per_kilogram"`, and `WorldOceanAtlasTools.convert_to_SI_unit!`
# # does `χ_3D .*= ustrip(upreferred(1.0u"μmol/kg"))`, which gives
# # `1.0e-6 mol/kg` — kg is already SI, so `upreferred` does **not**
# # convert to `mol/m³`. The returned array is therefore in `mol/kg`
# # (mean ≈ 1.6e-6), a factor of ~1000 smaller than the AIBECS model
# # state. We close the gap with a nominal seawater density before
# # anything downstream consumes these arrays.

# using WorldOceanAtlasTools
# μDIPobs3D, σ²DIPobs3D = WorldOceanAtlasTools.fit_to_grid(grd, "PO₄")
# iwet = iswet(grd)

# ρ_factor = ustrip(upreferred(1.035u"kg/L"))   # = 1035.0 kg/m³
# μDIPobs  = μDIPobs3D[iwet]  .* ρ_factor       # mol/kg → mol/m³
# σ²DIPobs = σ²DIPobs3D[iwet] .* ρ_factor^2     # (mol/kg)² → (mol/m³)²
# v = volumevec(grd)

# # ## Build the objective `(f, ∇ₓf)`

# # We weight the DIP misfit at 1, ignore POP (no direct observations),
# # and add a small parameter-prior penalty.

# ωs = (1.0, 0.0)
# ωp = 1.0e-4
# μx  = (μDIPobs,  missing)
# σ²x = (σ²DIPobs, missing)
# f, ∇ₓf = f_and_∇ₓf(ωs, μx, σ²x, v, ωp, OptParams)

# # ## F-1 memory cache + warm-started linear solve

# # `F1Method.initialize_mem` solves the steady-state problem once and
# # stashes the factorised Jacobian plus the state-parameter sensitivity
# # `∇s`; subsequent objective / gradient / Hessian calls reuse this
# # cache whenever the parameter vector has not moved. The
# # `F1MethodOptimizationExt` (loaded the moment `using Optimization`
# # runs) wraps the cached triple into a
# # `SciMLBase.OptimizationFunction`, exactly the way `LinearSolve.jl`
# # wraps UMFPACK into a SciML `LinearProblem` solver.
# #
# # We route F1Method's internal linear algebra through `LinearSolve.jl`
# # via the `linsolve` kwarg so the factorised Jacobian lives in a
# # `LinearSolve.LinearCache` we can share with AIBECS's inner Newton
# # solver. F1Method's `mem.linear_cache` carries a matrix RHS (the
# # parameter Jacobian `-∇ₚF`); AIBECS needs a vector RHS, so we keep a
# # sibling cache (`AIBECScache`) and share the UMFPACK factors via
# # `cacheval`. Threading this through *every* downstream `CTKAlg`
# # solve makes the cost-surface sweep below ~10× faster. When the
# # warm-started factors get too stale, CTKAlg's Shamanskii ratio check
# # + Armijo line search refresh them automatically.

# using F1Method, Optimization, OptimizationOptimJL, ADTypes
# using LinearSolve

# solver_kwargs = (; maxItNewton = 15)

# λ₀ = p2λ(p₀)
# mem = F1Method.initialize_mem(F, ∇ₓf, x0, λ₀, CTKAlg();
#     ad = AutoForwardDiff(),
#     linsolve = UMFPACKFactorization(),
#     solver_kwargs...)

# AIBECScache = init(
#     LinearProblem(F.jac(x0, λ2p(OptParams, λ₀)), similar(x0)),
#     UMFPACKFactorization(),
# )
# AIBECScache.cacheval = mem.linear_cache.cacheval
# AIBECScache.isfresh  = false

# solver_kwargs = (; solver_kwargs..., linear_cache = AIBECScache)

# optfn = F1Method.optimization_function(f, F, ∇ₓf, mem, CTKAlg(); solver_kwargs...)

# # ## Run the optimiser

# # `maxiters = 20` is plenty to watch the outer Newton trust-region
# # trace drive the gradient norm down by several decades. We attach a
# # callback that records the parameter vector at each iteration so the
# # diagnostics figure below can overlay the Newton trajectory on the
# # cost-surface contourf.

# λ_trajectory = Vector{Vector{Float64}}()
# push!(λ_trajectory, copy(λ₀))
# trajectory_cb = function (state, _)
#     push!(λ_trajectory, copy(state.u))
#     return false
# end

# prob = OptimizationProblem(optfn, copy(λ₀))
# results = solve(prob, NewtonTrustRegion();
#     maxiters = 20, callback = trajectory_cb, show_trace = true)
# p_opt = λ2p(OptParams, results.u)

# # ## Cost-function sweep around the trajectory

# # To visualise the loss landscape we evaluate `f(x_ss, λ)` on a 2D
# # grid in λ-space centred on the Newton trajectory, padded by ±0.2.
# # Each cell solves the steady-state PO₄ system independently from the
# # same cold-start `x0`, warm-started by the shared UMFPACK cache from
# # the previous section. The retcode is ignored — even when Newton
# # hits `maxItNewton` without fully converging, the residual is small
# # enough that `f(u, λ)` is a fine approximation of the true cost for a
# # coarse contour. With 10×10 = 100 solves this typically takes a
# # couple of minutes.

# n_grid = 10
# pad    = 0.2
# λ1_lo, λ1_hi = extrema(λ[1] for λ in λ_trajectory) .+ (-pad, pad)
# λ2_lo, λ2_hi = extrema(λ[2] for λ in λ_trajectory) .+ (-pad, pad)
# λ1_vec = range(λ1_lo, λ1_hi; length = n_grid)
# λ2_vec = range(λ2_lo, λ2_hi; length = n_grid)

# F_grid = Matrix{Float64}(undef, n_grid, n_grid)
# for j in 1:n_grid, i in 1:n_grid
#     λ = [λ1_vec[i], λ2_vec[j]]
#     sol = solve(SteadyStateProblem(F, x0, λ), CTKAlg(); solver_kwargs...)
#     F_grid[i, j] = f(sol.u, λ)
# end

# # ## Diagnostics — a single composed Makie figure
# #
# # We pull in `CairoMakie` (the headless backend used elsewhere in the
# # docs build — see the [Makie how-to](@ref plots-makie)) and bundle
# # every diagnostic into one `Figure` with four nested `GridLayout`s:
# #
# # - **G1** (top-left) — density-coloured scatter of optimised vs WOA
# #   PO₄: each wet box is plotted at `(WOA, optimised)`, coloured by a
# #   bivariate KDE of the point cloud (`cgrad(:viridis; rev = true)` so
# #   dense regions read dark), with a 1:1 reference line.
# # - **G2** (top-right) — `contourf!` of `√f(p)` in `(τ_DIP, τ_POP)`
# #   space on log axes, with the recorded Newton trajectory drawn on
# #   top (red line + markers, white-star initial point, red-star
# #   optimum).
# # - **G3** (middle row) — horizontal slices at 1000 m: optimised | WOA
# #   | difference, with `colorrange = (0, 3)` μM under viridis for the
# #   first two and `(-3, 3)` μM under `cgrad(:RdBu; rev = true)` for
# #   the difference.
# # - **G4** (bottom block) — 3-row × 4-column basin diagnostics grid:
# #   - Cols 1–3: per-basin (Pacific / Atlantic / Indian, each extended
# #     to the pole via its Southern Ocean sector) zonal-mean
# #     `contourf!`.
# #   - Col 4: a profile axis with one line per basin (depth profile
# #     across the basin's wet boxes).
# #   - Rows: optimised / WOA / optimised − WOA. The model and WOA rows
# #     share the `(0, 3)` μM viridis colorbar; the diff row gets its
# #     own `(-3, 3)` μM `:RdBu` (reversed) colorbar. Row identity is
# #     marked by a rotated `Label` in column 0.

using CairoMakie
using OceanBasins
using KernelDensity
CairoMakie.activate!(type = "png")

# Solve once more at the optimum to get the optimised DIP field.
prob_opt = SteadyStateProblem(F, x0, p_opt)
DIPopt, _ = state_to_tracers(solve(prob_opt, CTKAlg(); solver_kwargs...).u, grd)

# Strip units so the plotting calls work with plain numbers in μM.
DIPopt_μM = ustrip.(μM, DIPopt  .* (mol/m^3))
μWOA_μM   = ustrip.(μM, μDIPobs .* (mol/m^3))
diff_μM   = DIPopt_μM .- μWOA_μM

depth = ustrip.(grd.depth)
lat   = ustrip.(grd.lat)
lon   = ustrip.(grd.lon)
maxdepth_round = ceil(maximum(depth) / 1000) * 1000

# G2 — convert sweep grid and trajectory from λ-space to physical
# units. `λ2p` handles the inverse bijector per-parameter.
τDIP_grid = [λ2p(OptParams, [λ1, λ₀[2]]).τ_DIP for λ1 in λ1_vec]
τPOP_grid = [λ2p(OptParams, [λ₀[1], λ2]).τ_POP for λ2 in λ2_vec]
τDIP_traj = [λ2p(OptParams, λ).τ_DIP for λ in λ_trajectory]
τPOP_traj = [λ2p(OptParams, λ).τ_POP for λ in λ_trajectory]

# G3 — horizontal slices at 1000 m.
sl_opt  = AIBECS.horizontalslice(DIPopt_μM, grd; depth = 1000)
sl_obs  = AIBECS.horizontalslice(μWOA_μM,   grd; depth = 1000)
sl_diff = sl_opt .- sl_obs

# G4 — per-basin (extended-to-pole) masks for zonal means + profiles.
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

# Shared colour scales — fixed at user-requested ranges so every map /
# zonal panel uses the same bins.
dip_levels  = 0:0.25:3.0
dip_cmap    = cgrad(:viridis, length(dip_levels) - 1; categorical = true)
diff_levels = -3.0:0.5:3.0
diff_cmap   = cgrad(:RdBu, length(diff_levels) - 1; categorical = true, rev = true)

# Tick / label helpers.
latticks = (-90:30:90, ["90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"])
lonticks = (0:60:360, ["0°", "60°E", "120°E", "180°", "120°W", "60°W", "0°"])
gcdepthticks = (0:1000:maxdepth_round, string.(Int.(0:1000:maxdepth_round)))
plus_minus_ticks = (-3:1:3, ["−3", "−2", "−1", "0", "+1", "+2", "+3"])

fig = Figure(size = (1600, 1700))
g1 = fig[1, 1]   = GridLayout()
g2 = fig[1, 2]   = GridLayout()
g3 = fig[2, 1:2] = GridLayout()
g4 = fig[3, 1:2] = GridLayout()
rowsize!(fig.layout, 1, Auto(1.0))
rowsize!(fig.layout, 2, Auto(0.85))
rowsize!(fig.layout, 3, Auto(2.4))

## ── G1: density-coloured scatter ────────────────────────────────
ax_sc = Axis(g1[1, 1];
    xlabel = "WOA PO₄ (μM)", ylabel = "optimised PO₄ (μM)",
    xgridvisible = false, ygridvisible = false,
    title = "optimised vs WOA (KDE-coloured)")
## `kde` builds the gridded density; wrap it once in `InterpKDE` and
## broadcast on its underlying interpolator. The naïve
## `[pdf(d_kde, xi, yi) for (xi, yi) in ...]` rebuilds the InterpKDE
## (incl. the quadratic B-spline extrapolator) on *every* call —
## ~hundreds of seconds for ~84k points. Broadcasting `ik.itp` once
## is the same numerical result and runs in ~30 ms here.
d_kde = kde((μWOA_μM, DIPopt_μM))
ik = InterpKDE(d_kde)
c_dens = ik.itp.(μWOA_μM, DIPopt_μM)
sc = scatter!(ax_sc, μWOA_μM, DIPopt_μM;
    markersize = 2, color = c_dens,
    colormap = cgrad(:viridis; rev = true))
let xy_max = max(maximum(μWOA_μM), maximum(DIPopt_μM))
    lines!(ax_sc, [0, xy_max], [0, xy_max];
        color = :black, linestyle = :dash, label = "1:1")
end
axislegend(ax_sc; position = :lt, framevisible = false)
Colorbar(g1[2, 1], sc; vertical = false, flipaxis = false,
    label = "point density (KDE)")

## ── G2: cost surface + Newton trajectory ────────────────────────
ax_cs = Axis(g2[1, 1];
    xlabel = "τ_DIP (d)", ylabel = "τ_POP (d)",
    xscale = log10, yscale = log10,
    title = "√f(p) cost surface + Newton trajectory")
co_cs = contourf!(ax_cs, τDIP_grid, τPOP_grid, sqrt.(F_grid);
    levels = 20, colormap = :viridis)
lines!(ax_cs, τDIP_traj, τPOP_traj; color = :red, linewidth = 2)
scatter!(ax_cs, τDIP_traj, τPOP_traj;
    color = :red, marker = :circle, markersize = 7,
    strokecolor = :white, strokewidth = 0.5)
scatter!(ax_cs, [τDIP_traj[1]], [τPOP_traj[1]];
    marker = :star5, markersize = 22, color = :white,
    strokecolor = :black, strokewidth = 1.5, label = "p₀")
scatter!(ax_cs, [τDIP_traj[end]], [τPOP_traj[end]];
    marker = :star5, markersize = 22, color = :red,
    strokecolor = :black, strokewidth = 1.5, label = "p_opt")
axislegend(ax_cs; position = :rt, framevisible = false)
Colorbar(g2[1, 2], co_cs; label = "√f(p)")

## ── G3: horizontal slices at 1000 m ─────────────────────────────
function map_axis(grid_pos; title = "")
    Axis(grid_pos;
        backgroundcolor = :lightgray,
        xgridvisible = false, ygridvisible = false,
        xticks = lonticks, yticks = latticks,
        title = title,
        limits = (0, 360, -90, 90))
end

ax_sl_m = map_axis(g3[1, 1]; title = "optimised")
co_sl_m = contourf!(ax_sl_m, lon, lat, sl_opt';
    levels = dip_levels, colormap = dip_cmap, nan_color = :lightgray,
    extendlow = dip_cmap[1], extendhigh = dip_cmap[end])
ax_sl_o = map_axis(g3[1, 2]; title = "WOA")
contourf!(ax_sl_o, lon, lat, sl_obs';
    levels = dip_levels, colormap = dip_cmap, nan_color = :lightgray,
    extendlow = dip_cmap[1], extendhigh = dip_cmap[end])
ax_sl_d = map_axis(g3[1, 3]; title = "optimised − WOA")
co_sl_d = contourf!(ax_sl_d, lon, lat, sl_diff';
    levels = diff_levels, colormap = diff_cmap, nan_color = :lightgray,
    extendlow = diff_cmap[1], extendhigh = diff_cmap[end])

linkyaxes!(ax_sl_m, ax_sl_o, ax_sl_d)
hideydecorations!(ax_sl_o; ticklabels = true, label = false, ticks = false, grid = false)
hideydecorations!(ax_sl_d; ticklabels = true, label = false, ticks = false, grid = false)

cb_g3_v = Colorbar(g3[2, 1:2], co_sl_m; vertical = false, flipaxis = false,
    label = "PO₄ @ 1000 m (μM)", ticks = 0:1:3)
cb_g3_v.width = Relative(0.7)
cb_g3_d = Colorbar(g3[2, 3], co_sl_d; vertical = false, flipaxis = false,
    label = "optimised − WOA (μM)", ticks = plus_minus_ticks)
cb_g3_d.width = Relative(0.85)

## ── G4: per-basin zonal means + per-row basin profile (3×4) ─────
function zonal_axis_g4(grid_pos; ylabel = "", title = "")
    Axis(grid_pos;
        backgroundcolor = :lightgray,
        xgridvisible = false, ygridvisible = false,
        xticksmirrored = true, yticksmirrored = true,
        xticks = latticks, yticks = gcdepthticks,
        yreversed = true,
        title = title,
        limits = (-90, 90, 0, maxdepth_round),
        ylabel)
end

row_specs = [
    (label = "model",         field = DIPopt_μM, levels = dip_levels,  cmap = dip_cmap),
    (label = "WOA",           field = μWOA_μM,   levels = dip_levels,  cmap = dip_cmap),
    (label = "model − WOA",   field = diff_μM,   levels = diff_levels, cmap = diff_cmap),
]

## Representative contourf handles for the two colorbars below — we
## fill them inside the loop. Using `Ref{Any}` (instead of plain
## `co_g4_v = nothing`) because Documenter `@example` blocks run at
## module scope, where bare `co_g4_v = co` inside a `for` would shadow
## the outer binding; `Ref` is mutated via `[]`, which sidesteps the
## issue.
co_g4_v = Ref{Any}()
co_g4_d = Ref{Any}()

for (r, spec) in enumerate(row_specs)
    for (c, basin) in enumerate(basin_specs)
        ax = zonal_axis_g4(g4[r, c];
            ylabel = c == 1 ? "depth (m)" : "",
            title  = r == 1 ? basin.name : "")
        zm = AIBECS.zonalmean(spec.field, grd, basin.mask)
        co = contourf!(ax, lat, depth, zm;
            levels = spec.levels, colormap = spec.cmap, nan_color = :lightgray,
            extendlow = spec.cmap[1], extendhigh = spec.cmap[end])
        r == 1 && c == 1 && (co_g4_v[] = co)
        r == 3 && c == 1 && (co_g4_d[] = co)
        c > 1 && hideydecorations!(ax;
            ticklabels = true, label = true, ticks = false, grid = false)
        r < 3 && hidexdecorations!(ax;
            ticklabels = true, label = true, ticks = false, grid = false)
    end
    ## Profile axis in col 4 — one line per basin, depth on the right
    ## y-axis (the zonal column already shows depth on the left).
    ax_p = Axis(g4[r, 4];
        xgridvisible = false, ygridvisible = false,
        xlabel = r == 3 ? "PO₄ (μM)" : "",
        ylabel = "",
        yreversed = true,
        yticks = gcdepthticks,
        title = r == 1 ? "by basin" : "",
        limits = (nothing, (0, maxdepth_round)))
    for basin in basin_specs
        prof = AIBECS.horizontalmean(spec.field, grd, basin.mask)
        lines!(ax_p, prof, depth;
            color = basin.color, linewidth = 2, label = basin.name)
    end
    r == 1 && axislegend(ax_p; position = :rb, framevisible = false)
    r < 3 && hidexdecorations!(ax_p;
        ticklabels = true, label = true, ticks = false, grid = false)
end

## Row identity labels in column 0.
Label(g4[1, 0], "model";        rotation = π/2, font = :bold, fontsize = 16, tellheight = false)
Label(g4[2, 0], "WOA";          rotation = π/2, font = :bold, fontsize = 16, tellheight = false)
Label(g4[3, 0], "model − WOA";  rotation = π/2, font = :bold, fontsize = 16, tellheight = false)

## Shared vertical colorbars — viridis spans rows 1+2, RdBu sits under row 3.
cb_g4_v = Colorbar(g4[1:2, 5], co_g4_v[];
    label = "PO₄ (μM)", ticks = 0:1:3)
cb_g4_v.height = Relative(0.85)
cb_g4_d = Colorbar(g4[3, 5], co_g4_d[];
    label = "model − WOA (μM)", ticks = plus_minus_ticks)
cb_g4_d.height = Relative(0.85)

## Column sizing — narrow row-label / colorbar cols, equal zonal cols,
## slightly narrower profile col.
colsize!(g4, 0, Auto(0.15))
colsize!(g4, 1, Auto(2.0))
colsize!(g4, 2, Auto(2.0))
colsize!(g4, 3, Auto(2.0))
colsize!(g4, 4, Auto(1.2))
colsize!(g4, 5, Auto(0.2))

## Breathing room between the four outer grid blocks.
colgap!(fig.layout, 1, 30)
rowgap!(fig.layout, 1, 30)
rowgap!(fig.layout, 2, 30)

fig

# Side-by-side parameter tables, initial vs. optimised:

AIBECS.table(p₀)

#-

AIBECS.table(p_opt)

# ## Where to next

# This example only fits two parameters of a single-tracer mismatch
# objective on a moderate-resolution circulation. AIBECS also exposes
# the table-form `f_and_∇ₓf(ωs, ωp, grd, obs, ::Type{P})` for fitting
# directly against `WorldOceanAtlasTools.observations(...)` tables and
# the `f_and_∇ₓf(ωs, ωp, grd, modify, obs, ::Type{P})` form for
# tracer-derived quantities (isotopes via δ-conversion, etc.); the
# multi-tracer, multi-observation, multi-iteration fit lives in
# [GNOM](https://github.com/MTEL-USC/GNOM).

# Related, but not used here:
#
# - [`SciMLSensitivity.SteadyStateAdjoint`](https://docs.sciml.ai/SciMLSensitivity/stable/)
#   computes a per-call adjoint gradient, with no cross-call caching and
#   no Hessian — appropriate when the steady-state solve dominates and
#   you only need the gradient.
# - [`ImplicitDifferentiation.jl`](https://github.com/gdalle/ImplicitDifferentiation.jl)
#   provides general implicit-function differentiation; no F-1 Hessian
#   today.
#
# Both are good upstream PR targets to slot underneath the AIBECS
# objective / gradient surface in the future.

# ## Caveats and possible future directions

# The setup above is fairly bespoke: AIBECS' parameter struct + a
# hand-written `prior` method + F1Method's cached adjoint + a sibling
# `LinearSolve` cache + the `Optimization.jl` frontend. It works, but
# it is not the only — nor necessarily the best — way to attack
# parameter estimation for steady-state biogeochemical models. A few
# directions where AIBECS could grow:
#
# - **Probabilistic programming via [Turing.jl](https://turinglang.org/).**
#   Today we declare priors with `@prior` (or by overloading `prior`
#   directly) and the optimiser minimises `-log posterior` up to a
#   constant. Wiring AIBECS parameter structs into a Turing `@model`
#   would unlock full Bayesian inference (MCMC, variational families)
#   over the steady-state solution, with the same priors driving both
#   regimes. The F-1 gradient/Hessian would compose naturally with
#   Turing's HMC / NUTS samplers.
# - **Adopt the SciML parameter-estimation idioms.** The SciML stack
#   (`SciMLSensitivity`, `DiffEqParamEstim`, `Optimization`) has a
#   well-trodden pipeline for fitting `ODEProblem` / `SteadyStateProblem`
#   parameters from data. Right now AIBECS bolts F1Method onto
#   `Optimization` by hand. Restructuring the AIBECS objective surface
#   so it slots into SciML's `build_loss_objective`-style helpers would
#   let users mix AIBECS models with the rest of the SciML param-fitting
#   ecosystem without bespoke glue.
# - **First-class observation tables.** GNOM's table-form
#   `f_and_∇ₓf(ωs, ωp, grd, obs, ::Type{P})` is the right shape for
#   GEOTRACES-style cruise data; promoting it (and adding more
#   adapters) would replace the volume-mean / variance shortcut used
#   here with proper interpolation-aware misfits.
#
# Each of these would be a meaningful upstream PR; the current docs
# example is closer to a worked sketch than a reference workflow.
