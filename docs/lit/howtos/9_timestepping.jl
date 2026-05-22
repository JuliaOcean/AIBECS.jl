#---------------------------------------------------------
# # [Time-stepping](@id timestepping)
#---------------------------------------------------------

# AIBECS is built around the steady-state form of the tracer equation
#
# $$\partial_t \boldsymbol{x} = \boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = \boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p}) - \mathbf{T}(\boldsymbol{p}) \, \boldsymbol{x} = 0,$$
#
# and the default solver, `CTKAlg()`, finds the steady state directly via
# Newton iteration. For some workflows you instead want to *march in time* —
# producing a trajectory from an initial condition to a finite end time, or
# spinning a model up through transient damping rather than a single Newton
# solve. AIBECS supports this by handing the same [`AIBECSFunction`](@ref)
# you already build for the steady-state problem to any SciML time-stepping
# integrator.

# Two helpers do the glue:
#
# - [`AIBECS.odeproblem`](@ref) wraps an `AIBECSFunction` into a
#   `SciMLBase.ODEProblem` with the sparse Jacobian buffer attached
#   (otherwise stiff solvers blow up the memory with a dense buffer).
# - [`AIBECS.splitodeproblem`](@ref) wraps an [`AIBECSSplitFunction`](@ref)
#   into a `SciMLBase.SplitODEProblem` so IMEX schemes can treat the linear
#   transport implicitly and the nonlinear sources explicitly.

# The actual time steppers come from
# [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) and its
# sub-packages. AIBECS does not re-export them — load the sub-package that
# carries the scheme you want:
#
# - `OrdinaryDiffEqTsit5` — `Tsit5` (explicit 5(4) RK)
# - `OrdinaryDiffEqSDIRK` — `ImplicitEuler`, `Trapezoid`, `TRBDF2`, `KenCarp4`
# - `OrdinaryDiffEqBDF` — `SBDF2`, `FBDF`, `IMEXEuler`
# - `OrdinaryDiffEqIMEXMultistep` — `CNAB2`, `CNLF2`

using AIBECS
using OrdinaryDiffEqSDIRK   # ImplicitEuler
using OrdinaryDiffEqBDF     # SBDF2
using LinearSolve           # UMFPACKFactorization
using Unitful: yr, Myr, kyr

# We use the toy `Primeau_2x2x2` circulation so the docs build stays fast;
# any AIBECS circulation works the same way.

grd, T = Primeau_2x2x2.load()
nb = sum(iswet(grd))
z = depthvec(grd)

# ## A model to play with: ideal age

# The ideal-age tracer satisfies
#
# $$\partial_t \boldsymbol{a} + \mathbf{T} \, \boldsymbol{a} = 1 - \frac{\boldsymbol{a}}{\tau} \, (\boldsymbol{z} \le z_0),$$
#
# i.e. it ages at 1 s/s in the interior and is restored to zero near the
# surface. This is **fully linear in `a`** — which matters below for the
# IMEX path. `G` and the linear part `L` are the same; the nonlinear part
# `NL` is identically zero.

function G(x, p)
    @unpack τ, z₀ = p
    return @. 1 - x / τ * (z ≤ z₀)
end
L_lin(x, p) = G(x, p)
NL_zero(x, p) = zero(x)

struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end
p = IdealAgeParameters(1.0, z[1])

fun = AIBECSFunction(T, G)
splitfun = AIBECSSplitFunction(T, L_lin, NL_zero, nb)

# ## Reference: steady state via `CTKAlg`

# This is what we'll cross-check the time integrations against.

x₀ = zeros(nb)
sol_ss = solve(SteadyStateProblem(fun, x₀, p), CTKAlg())
nothing # hide

# ## 1. Explicit time-stepping (toy use)

# Any explicit scheme can integrate `AIBECSFunction` directly via
# `ODEProblem` — no AIBECS helper needed. This works because
# `AIBECSFunction` already implements the `f(u, p, t)` and `jac(u, p, t)`
# signatures SciML expects.

# !!! warning "CFL on real circulations"
#     Explicit schemes are CFL-limited by the transport operator. On a
#     production grid like OCIM2 the stable `dt` is on the order of
#     **hours**, so anything beyond a brief diagnostic transient becomes
#     prohibitively many steps. Use implicit or IMEX for spin-up.

using OrdinaryDiffEqTsit5
t_end_short = ustrip(upreferred(1.0yr))
prob_explicit = ODEProblem(fun, x₀, (0.0, t_end_short), p)
sol_explicit = solve(prob_explicit, Tsit5())
nothing # hide

# ## 2. Fully implicit time-stepping

# For spin-up, the workhorse is `ImplicitEuler` (L-stable, robust, first
# order) — useful when you don't care about accuracy along the trajectory,
# just about reaching the steady state. For trajectory accuracy use
# `Trapezoid` (Crank–Nicolson, 2nd order) or higher-order SDIRK methods
# like `KenCarp4` / `TRBDF2`.

# Wrap with `AIBECS.odeproblem` so the sparse Jacobian prototype is
# attached — otherwise OrdinaryDiffEq's stiff path allocates a dense `N×N`
# Jacobian buffer (catastrophic at production scale).

t_end = ustrip(upreferred(10.0Myr))
dt = ustrip(upreferred(100.0kyr))

prob_implicit = AIBECS.odeproblem(fun, x₀, (0.0, t_end), p)
sol_implicit = solve(
    prob_implicit, ImplicitEuler(linsolve = UMFPACKFactorization());
    dt = dt, saveat = [t_end],
)
maximum(abs, sol_implicit.u[end] - sol_ss.u) / maximum(abs, sol_ss.u)

# ### Constant-Jacobian models: factorization reuse

# Ideal age has a Jacobian `∂F/∂x = -L_diag - T` that does not depend on
# `x`. Pass `jac_constant_in_u = true` so the helper caches the Jacobian
# once and signals to the integrator that the W matrix can be factorised
# once and reused across all steps. This is the single biggest win the
# helpers expose, and it applies to any model where `G` is affine in `x`
# (ideal age, radiocarbon, linear restoring, …).

prob_constjac = AIBECS.odeproblem(fun, x₀, (0.0, t_end), p; jac_constant_in_u = true)
sol_constjac = solve(
    prob_constjac, ImplicitEuler(linsolve = UMFPACKFactorization());
    dt = dt, saveat = [t_end],
)
maximum(abs, sol_constjac.u[end] - sol_ss.u) / maximum(abs, sol_ss.u)

# !!! warning "Only set `jac_constant_in_u = true` when the math actually permits"
#     This is a user assertion, not a fact we verify. If `G` depends on `x`
#     non-affinely (e.g. Michaelis–Menten uptake), passing `true` will give
#     the integrator a stale W matrix and the solution will be silently
#     wrong. Default `false` is the safe choice.

# ## 3. IMEX time-stepping

# When transport (always linear) plus part of the sources-and-sinks are
# linear-in-`u`, IMEX schemes can split the equation so the *stiff linear*
# half is treated implicitly (one linear solve per step) and the *non-stiff
# nonlinear* half is treated explicitly. The win: each step is dominated by
# one back-solve against a factorisation that can be **reused** across the
# whole integration as long as `dt` is constant.

# Build the split with [`AIBECSSplitFunction`](@ref) — it mirrors
# [`AIBECSFunction`](@ref) but takes separate tuples for the linear sources
# `Ls` and nonlinear sources `NLs`. Wrap with `AIBECS.splitodeproblem` to
# attach the (constant in `u`) Jacobian of the implicit half:

sprob = AIBECS.splitodeproblem(splitfun, x₀, (0.0, t_end), p)
sol_imex = solve(
    sprob, SBDF2(linsolve = UMFPACKFactorization());
    dt = dt, saveat = [t_end],
)
maximum(abs, sol_imex.u[end] - sol_ss.u) / maximum(abs, sol_ss.u)

# !!! warning "IMEX requires all stiff dynamics to live in `L`"
#     IMEX schemes treat `f₂ = NL` *explicitly*, so any stiff timescale in
#     `NL` imposes a CFL-like bound `dt ≲ τ_stiff`. If your nonlinear sources
#     include stiff reactions (e.g. fast Michaelis–Menten uptake with
#     `τ ≈ 30 days`), an IMEX run with `dt = 1 kyr` will diverge silently —
#     the integrator returns `Success` but the answer is nonsense. For
#     models with nonlinear stiff sources, prefer fully implicit methods
#     (`ImplicitEuler`, `KenCarp4`) over IMEX. IMEX shines when **all**
#     stiff dynamics is linear: ideal age, radiocarbon, linear restoring,
#     linear remineralisation chains.

# ## Choosing `dt`

# Heuristic: start with `dt ≈ τ_slowest / 10`, where `τ_slowest` is the
# longest relaxation timescale you actually care about resolving. For
# tracer transport in the global ocean `τ_slowest` is typically the
# overturning timescale (~10³ years), so `dt ≈ 100 yr` is a reasonable
# starting point.
#
# - Smaller `dt` ⇒ better trajectory accuracy, more steps.
# - Larger `dt` ⇒ faster spin-up but loses time-accuracy. With
#   `ImplicitEuler` this is fine — it just damps to steady state.
# - For multistep IMEX (`SBDF2`, `CNAB2`) keep `dt` *constant* —
#   adaptive `dt` forces a refactorisation at every change, defeating the
#   purpose. For adaptive stiff stepping use `KenCarp4` or `TRBDF2`.

# ## Cross-check against `CTKAlg`

# A useful sanity test for any time-stepped result: re-solve via
# `CTKAlg` and confirm the trajectory's endpoint matches the steady state
# within the integrator's accuracy tolerance.

@assert maximum(abs, sol_implicit.u[end] - sol_ss.u) ≤ 1e-2 * maximum(abs, sol_ss.u)
@assert maximum(abs, sol_imex.u[end]    - sol_ss.u) ≤ 1e-2 * maximum(abs, sol_ss.u)
