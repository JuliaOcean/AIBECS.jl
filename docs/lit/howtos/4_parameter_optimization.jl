#---------------------------------------------------------
# # [Parameter optimisation](@id parameter-optimization)
#---------------------------------------------------------

# In this how-to we fit three parameters of the coupled PO₄–POP model
# (`τ_DIP`, `τ_POP`, `k`) to PO₄ observations from
# [WorldOceanAtlasTools.jl](https://github.com/briochemc/WorldOceanAtlasTools.jl)
# on the OCCA circulation. We use the SciML [Optimization.jl][1] frontend,
# with [F1Method.jl][2] supplying the cached steady-state gradient and
# Hessian.
#
# The API surface we touch:
#
# - `@flattenable @bounds @logscaled @prior` — declare which parameters
#   are optimised and supply their priors;
# - [`p2λ`](@ref) / [`λ2p`](@ref) — transform a parameter struct to/from
#   an unconstrained flat vector;
# - [`f_and_∇ₓf`](@ref) — build the volume-weighted mismatch objective
#   and its analytical state Jacobian;
# - `F1Method` — cached steady-state solve plus adjoint gradient and
#   Hessian ([Pasquier et al., 2019](https://doi.org/10.5194/gmd-12-3537-2019));
# - `Optimization.jl` — the SciML optimisation frontend, in the same way
#   `LinearSolve.jl` is the SciML linear-solver frontend used elsewhere
#   in AIBECS.
#
# See the [coupled PO₄–POP tutorial](@ref P-model) for the model
# equations and [GNOM](https://github.com/MTEL-USC/GNOM) for the
# converged multi-tracer / isotope fit this example is a thin slice of.
#
# [1]: https://docs.sciml.ai/Optimization/stable/
# [2]: https://github.com/briochemc/F1Method.jl

# ## Re-build the P-model on OCCA

# We reuse the source/sink terms from the P-model tutorial verbatim, just
# swapping the circulation to OCCA so the resolution stays moderate
# (~170k wet boxes × 2 tracers).

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
# (canonical example at `src/Parameters.jl`); the columns therefore are
# `initial_value | units | flattenable | bounds | logscaled | prior`.

import AIBECS: @initial_value, initial_value
import AIBECS: @units, units
import AIBECS: @flattenable, flattenable
import AIBECS: @bounds, bounds
import AIBECS: @logscaled, logscaled
import AIBECS: prior
using Unitful: m, d, s, yr, Myr, mol, mmol, μmol, μM
using Distributions
using Bijectors    # activates the `AIBECSBijectorsExt` extension (p2λ, λ2p)
using DataFrames   # activates the `AIBECSDataFramesExt` extension (AIBECS.table)

const ∞ = Inf

@initial_value @units @flattenable @bounds @logscaled struct OptParams{U} <: AbstractParameters{U}
    w₀::U      |   0.64 | m/d      | false | (  0,     ∞) | true
    w′::U      |   0.13 | m/d/m    | false | (  0,     ∞) | true
    τ_DIP::U   | 230.0  | d        | true  | (  1, 500.0) | true
    k::U       |   6.62 | μmol/m^3 | true  | (0.1,  10.0) | true
    z₀::U      |  80.0  | m        | false | (  0,     ∞) | false
    τ_POP::U   |   5.0  | d        | true  | (  0,     ∞) | true
    τ_geo::U   |   1.0  | Myr      | false | (  0,     ∞) | true
    DIP_geo::U |   2.12 | mmol/m^3 | false | (  0,     ∞) | false
end

# Derive priors from `(lb, ub) = bounds(T, s)` and the initial value, matching
# the convention in `test/parameters.jl` and GNOM: positive half-line gets a
# LogNormal, the real line a wide Normal, and finite-interval parameters get a
# LogitNormal stretched onto `(lb, ub)`.

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

# Three fields (`τ_DIP`, `k`, `τ_POP`) carry `flattenable = true`; the
# rest are held fixed at their initial value. We can render the full
# metadata table with [`AIBECS.table`](@ref):

AIBECS.table(p₀)

# ## State function in parameter-type form

# Using the `AIBECSFunction(Ts, Gs, nb, ::Type{P})` form lets the SciML
# solvers feed in a flat unconstrained vector `λ`; AIBECS will call
# `λ2p(OptParams, λ)` internally to reconstruct the struct.

nb = sum(iswet(grd))
F = AIBECSFunction((T_DIP, T_POP), (G_DIP, G_POP), nb, OptParams)
x0 = p₀.DIP_geo * ones(2nb)

# ## Load WOA PO₄ observations

# `WorldOceanAtlasTools.fit_to_grid` returns the mean and variance of
# WOA PO₄ regridded onto our `grd` as 3D fields; we restrict to wet
# boxes and grab the box volumes for the volume-weighted misfit.

using WorldOceanAtlasTools
μDIPobs3D, σ²DIPobs3D = WorldOceanAtlasTools.fit_to_grid(grd, "PO₄")
iwet = iswet(grd)
μDIPobs  = μDIPobs3D[iwet]
σ²DIPobs = σ²DIPobs3D[iwet]
v = volumevec(grd)

# ## Build the objective `(f, ∇ₓf)`

# We weight the DIP misfit at 1, ignore POP (no direct observations),
# and add a small parameter-prior penalty.

ωs = (1.0, 0.0)
ωp = 1.0e-4
μx  = (μDIPobs,  missing)
σ²x = (σ²DIPobs, missing)
f, ∇ₓf = f_and_∇ₓf(ωs, μx, σ²x, v, ωp, OptParams)

# ## F-1 memory cache and `OptimizationFunction`

# `F1Method.initialize_mem` solves the steady-state problem once and
# stashes the factorised Jacobian plus the state-parameter sensitivity
# `∇s`; subsequent objective / gradient / Hessian calls reuse this cache
# whenever the parameter vector has not moved. The
# `F1MethodOptimizationExt` (loaded the moment `using Optimization` runs)
# wraps the cached triple into a `SciMLBase.OptimizationFunction`,
# exactly the way `LinearSolve.jl` wraps UMFPACK into a SciML
# `LinearProblem` solver.

using F1Method, Optimization, OptimizationOptimJL, ADTypes

# CTKAlg's residual-norm tolerance and Newton iteration cap. Tightening
# `abstol` past the default makes the cached gradient/Hessian more
# accurate at the cost of inner Newton work per outer optimisation step.
# Both `initialize_mem` and `optimization_function` forward these kwargs
# straight to `solve(::SteadyStateProblem, ::CTKAlg; ...)`.

solver_kwargs = (; abstol = 1.0e-10, maxItNewton = 50)

λ₀ = p2λ(p₀)
mem = F1Method.initialize_mem(F, ∇ₓf, x0, λ₀, CTKAlg();
                              ad = AutoForwardDiff(), solver_kwargs...)
optfn = F1Method.optimization_function(f, F, ∇ₓf, mem, CTKAlg(); solver_kwargs...)

# ## Run the optimiser

# `maxiters = 20` is plenty to watch the outer Newton trust-region
# trace drive the gradient norm down by several decades. This is still
# a workflow demonstration rather than a converged fit — see
# [GNOM](https://github.com/MTEL-USC/GNOM) for the multi-tracer /
# multi-observation production version.

prob = OptimizationProblem(optfn, λ₀)
results = solve(prob, NewtonTrustRegion(); maxiters = 20, show_trace = true)
p_opt = λ2p(OptParams, results.u)

# ## Diagnostics

# Solve once more at the optimum, then compare the optimised DIP field
# against the WOA mean.

using Plots
prob_opt = SteadyStateProblem(F, x0, p_opt)
DIPopt, _ = state_to_tracers(solve(prob_opt, CTKAlg()).u, grd)
plothorizontalmean(DIPopt * (mol/m^3) .|> μM, grd, label = "optimized")
plothorizontalmean!(μDIPobs * (mol/m^3) .|> μM, grd, label = "WOA")

# Side-by-side parameter tables, initial vs. optimised:

AIBECS.table(p₀)

#-

AIBECS.table(p_opt)

# ## Where to next

# This example only fits three parameters of a single-tracer mismatch
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
# hand-written `prior` method + F1Method's cached adjoint + the
# `Optimization.jl` frontend. It works, but it is not the only — nor
# necessarily the best — way to attack parameter estimation for
# steady-state biogeochemical models. A few directions where AIBECS
# could grow:
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
