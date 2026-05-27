using ForwardDiff: Dual
using Unitful: s, Myr

"""
    CTKAlg{L} <: SciMLBase.AbstractSteadyStateAlgorithm

Marker algorithm that selects AIBECS's customised Newton–Chord–Shamanskii
solver when passed to `SciMLBase.solve(::SteadyStateProblem, ::CTKAlg)`.
The `linsolve` field is forwarded to `LinearSolve.LinearProblem`, so any
`SciMLLinearSolveAlgorithm` (UMFPACK, KLU, MKL Pardiso, …) can be plugged
in to swap the inner linear solver while keeping the Shamanskii nonlinear
logic.

```julia
prob = SteadyStateProblem(fun, x₀, p)
sol  = solve(prob, CTKAlg(); maxItNewton = 50)
sol  = solve(prob, CTKAlg(linsolve = KLUFactorization()))
```

The underlying engine is `AIBECS.NewtonChordShamanskii` (with Armijo line
search and lazy Jacobian refresh, after Kelley 2003).
"""
struct CTKAlg{L} <: SciMLBase.AbstractSteadyStateAlgorithm
    linsolve::L
end
CTKAlg(; linsolve = UMFPACKFactorization()) = CTKAlg(linsolve)

# Default `reltol` for CTKAlg, expressed as the inverse of a physical
# "steady-state" timescale τstop = 1 Myr in SI seconds. Since AIBECS state
# variables `x` carry SI units of (concentration | age | …) and the residual
# `F = ∂ₜx` carries (… / second), the ratio `‖F‖/‖x‖` is a rate in 1/s, and
# the criterion `‖F‖/‖x‖ ≤ reltol` (i.e. `‖x‖/‖F‖ ≥ τstop`) reads as "the
# system's typical change timescale exceeds τstop." 1 Myr is comfortably
# beyond any ocean-tracer transient yet still achievable at Float64 precision
# for typical AIBECS Jacobian conditioning.
const τstop_default = ustrip(s, 1Myr)        # 1 Myr in seconds ≈ 3.16e13
const reltol_default = 1 / τstop_default      # ≈ 3.17e-14

"""
    default_termination_condition() -> NamedTuple

Return AIBECS's recommended termination configuration as a `NamedTuple`
suitable for splatting into any `solve` call:

```julia
solve(prob, alg; AIBECS.default_termination_condition()..., maxiters = 50)
```

Produces the same stopping rule `CTKAlg` uses by default:
`NormTerminationMode` on the ∞-norm of the residual,
`reltol = 1 / τstop_default ≈ 3.17e-14` (i.e. steady-state timescale
≥ 1 Myr), and `abstol = 0`.

Splat it through to every solver that accepts `termination_condition` —
`CTKAlg` and the NonlinearSolve algorithms (`NewtonRaphson`,
`TrustRegion`, …) — so cross-solver comparisons share the same
convergence criterion. `NLsolveJL` is the exception: its wrapper does
not accept `termination_condition` and reports an L2 norm internally.
"""
function default_termination_condition()
    return (;
        termination_condition = NonlinearSolveBase.NormTerminationMode(
            Base.Fix1(maximum, abs)
        ),
        abstol = 0.0,
        reltol = reltol_default,
    )
end

"""
    solve(prob::SciMLBase.AbstractSteadyStateProblem,
          alg::CTKAlg;
          termination_condition = NormTerminationMode(Base.Fix1(maximum, abs)),
          abstol = 0.0,
          reltol = 1 / (1 Myr in seconds) ≈ 3.17e-14,
          preprint = "",
          maxItNewton = 50,
          linear_cache = nothing)

Solves `prob` using an AIBECS-customized Newton–Chord–Shamanskii method with
Armijo line search. Convergence is delegated to NonlinearSolveBase termination
modes.

The default is `NormTerminationMode` on the ∞-norm of the residual, which
terminates when either

    ‖F‖∞ ≤ abstol           (absolute floor)
    ‖F‖∞ ≤ reltol · ‖F + x‖∞  (relative; ≈ ‖F‖∞ ≤ reltol · ‖x‖∞ near steady state)

is satisfied. With the default `abstol = 0` the active branch is the relative
one. Why is this physically meaningful? Because`‖x‖/‖F‖` is the typical timescale
on which `x` is still changing, so `‖F‖/‖x‖ ≤ reltol` means that timescale
exceeds `τstop = 1 / reltol` seconds. The default `reltol ≈ 3.17e-14` corresponds
to `τstop = 1 Myr`, which is comfortably beyond any ocean-tracer transient.

Pass `abstol`/`reltol` to tighten or loosen, or pass any
`AbstractNonlinearTerminationMode` (e.g. `AbsNormSafeBestTerminationMode(norm)`,
`RelNormTerminationMode(norm)`) to swap criteria entirely.

A warning (`@warn`) is emitted if the solver returns a non-success `retcode`.

# Warm-starting via `linear_cache`

`linear_cache` (advanced) accepts a `LinearSolve.LinearCache` from a previous
solve. When supplied it replaces the internal `init(LinearProblem(…))` call,
so the caller's factorisation is reused on the first Newton iteration. This
matters in two settings:

1. Outer optimisers (e.g. `F1Method`, `Optimization.jl`) that already hold a
   factorised Jacobian for the current parameter point — feeding it back in
   skips the most expensive operation of the inner Newton solve.
2. Parameter sweeps where adjacent grid cells share most of the Jacobian's
   numerical content — the previous cell's factor is a useful warm-start.

The caller is responsible for setting `linear_cache.A`,
`linear_cache.cacheval`, and `linear_cache.isfresh` to a coherent triple
before the call. The cache's RHS must be a vector (length `nb`). Shamanskii's
ratio check + Armijo's line-search are the safety net if the supplied
factors turn out to be too stale — a bad warm-start slows convergence by at
most one Newton iteration before the factorisation is refreshed.
"""
function SciMLBase.solve(
        prob::SciMLBase.AbstractSteadyStateProblem,
        alg::CTKAlg;
        termination_condition = NonlinearSolveBase.NormTerminationMode(
            Base.Fix1(maximum, abs)
        ),
        abstol = 0.0,
        reltol = reltol_default,
        preprint = "",
        maxItNewton = 50,
        linear_cache = nothing,
    )
    # Define the functions according to SciMLBase.SteadyStateProblem type
    p = prob.p
    x0 = copy(prob.u0)
    dx = copy(x0)
    function F(x)
        return if SciMLBase.isinplace(prob)
            prob.f(dx, x, p)
            dx
        else
            prob.f(x, p)
        end
    end
    ∇ₓF(x) = prob.f.jac(x, p)
    F0 = F(x0)
    # Build the NonlinearSolveBase termination cache once. `prob` is just a
    # dispatch argument; SteadyStateProblem <: AbstractNonlinearProblem, so
    # this works without converting to a NonlinearProblem.
    tc_cache = NonlinearSolveBase.init(
        prob, termination_condition, F0, x0; abstol, reltol
    )
    # Compute `u_steady` and `resid` as per SciMLBase using my algorithm.
    # `alg.linsolve` is set eagerly to `UMFPACKFactorization()` by the default
    # constructor — see `CTKAlg(; linsolve = UMFPACKFactorization())` above —
    # so no `=== nothing` branch is needed here.
    x_steady, stats, retcode = NewtonChordShamanskii(
        F, ∇ₓF, x0, alg.linsolve, tc_cache;
        preprint = preprint, maxItNewton = maxItNewton, linear_cache = linear_cache
    )
    resid = F(x_steady)
    if !SciMLBase.successful_retcode(retcode)
        residual_inf_norm = maximum(abs, resid)
        @warn "CTKAlg did not converge to a successful retcode." retcode residual_inf_norm abstol reltol maxItNewton
    end
    # Map CTKAlg's internal counts onto `SciMLBase.NLStats` so the harness
    # reads stats uniformly across solvers. `nf` left at 0: F-call counting
    # would need a wrapper around `F` (and a separate one inside
    # `searchLineArmijo!`); not worth the complexity for now. CTKAlg does
    # exactly one linear solve per Newton step, so `nsolve = nsteps`.
    nl_stats = SciMLBase.NLStats(
        0,                          # nf — not tracked
        stats.jacobian_refreshes,   # njacs
        stats.jacobian_refreshes,   # nfactors — every Jacobian refresh triggers a refactor
        stats.iterations,           # nsolve — 1 linear solve per Newton step
        stats.iterations,           # nsteps
    )
    return SciMLBase.build_solution(
        prob, alg, x_steady, resid; retcode = retcode, stats = nl_stats,
    )
end

export solve, SteadyStateProblem, CTKAlg, default_termination_condition

# ── AD-through-the-solve via the implicit function theorem ───────────────────
# Mirror of `NonlinearSolveBase.nonlinearsolve_forwarddiff_solve` (defined in
# `NonlinearSolveBaseForwardDiffExt`) for AIBECS's `CTKAlg` route. SciML's IFT
# dispatch only fires for `NonlinearProblem`/`NonlinearLeastSquaresProblem`;
# `SteadyStateProblem` with a Dual-typed `p` does not hit it, so we add the
# equivalent dispatch here.
const DualSteadyStateProblem = SciMLBase.SteadyStateProblem{
    <:AbstractArray, iip, <:AbstractArray{<:Dual{T, V, P}},
} where {iip, T, V, P}

"""
    solve(prob::SteadyStateProblem{...,Dual,...}, alg::CTKAlg; ...)

When `prob.p` carries `ForwardDiff.Dual` numbers (i.e. AIBECS's `solve` is
invoked *inside* `ForwardDiff.gradient` / `ForwardDiff.hessian`), strip the
duals, run the primal solve via the existing CTKAlg method, and reconstruct
the dual `u*` via the implicit function theorem:

    ∂u*/∂p  =  -(∂f/∂u)⁻¹ ∂f/∂p     (at u*, p)

Reuses `prob.f.jac` for `∂f/∂u` (analytical/sparse from `AIBECSFunction`).
For nested duals (Hessian path), this dispatch fires recursively — one Dual
layer is peeled per recursion level, with the innermost solve hitting the
Float64 path.

!!! note "Production-scale Hessians"
    `ForwardDiff.hessian` over this dispatch costs `O(N_partials)` Newton
    solves and uses LinearSolve's `DualLinearCache` splitting (sparse all the
    way down). It is fine for small/medium circulations but not competitive
    with `F1Method`'s analytical gradient/Hessian formula, which reuses a
    cached steady state. Prefer `F1Method` for production-scale Hessians.

!!! warning "Partial-derivative accuracy"
    The primal Newton iteration terminates when `‖F(u)‖` meets `abstol`/`reltol`
    — only the *primal* `u*` is guaranteed to that tolerance. The IFT-reconstructed
    `∂u*/∂p` is evaluated at the not-quite-converged `u_steady` and inherits the
    primal residual: a sloppy primal produces sloppy partials, and the solver
    will not run extra iterations on the partials' behalf. To get derivatives to
    near-machine precision you must force the solver past its own stopping rule,
    e.g. `abstol = 0, reltol = 0` and let it grind out to `maxItNewton`. AIBECS's
    `test/ad_through_solve.jl` and the cross-check in F1Method's
    `bp-claude/ift-cross-check` branch show this pattern on toy problems.
    **This is not a practical strategy at production grid sizes** — each forced
    iteration costs a sparse factorisation. Use `F1Method`'s analytical
    gradient/Hessian for production-scale derivatives.

!!! warning "NonlinearSolveBase coupling"
    This dispatch depends on `NonlinearSolveBase.nodual_value`,
    `nonlinearsolve_∂f_∂p`, and `nonlinearsolve_dual_solution`. They are
    declared public via `@compat(public, …)` in NonlinearSolveBase but live
    in an internal-feeling part of the SciML stack. If they move or change
    signature, this dispatch breaks; a runtime `@warn` (below, gated to
    one-shot via `maxlog=1`) surfaces the dependency in user logs.
"""
function SciMLBase.solve(prob::DualSteadyStateProblem, alg::CTKAlg, args...; kwargs...)
    @warn """
        AIBECS's `solve(::SteadyStateProblem{…,Dual,…}, ::CTKAlg)` IFT dispatch \
        is exercising NonlinearSolveBase internals (`nodual_value`, \
        `nonlinearsolve_∂f_∂p`, `nonlinearsolve_dual_solution`). They are \
        declared public but are not on a stable SciML-wide API contract; \
        a breaking change upstream would break this dispatch. For \
        production-scale Hessians prefer `F1Method`.
        """ maxlog = 1

    # 1. Peel one Dual layer and remake the SteadyStateProblem at the value
    #    eltype.
    p_val = NonlinearSolveBase.nodual_value(prob.p)
    u0_val = NonlinearSolveBase.nodual_value(prob.u0)
    primal_prob = SciMLBase.remake(prob; p = p_val, u0 = u0_val)
    sol = solve(primal_prob, alg, args...; kwargs...)
    u_steady = sol.u

    # 2. IFT system: Jᵤ = ∂f/∂u (sparse, analytical via AIBECSFunction);
    #    Jₚ = ∂f/∂p (dense, ForwardDiff via NonlinearSolveBase helper).
    Jᵤ = prob.f.jac(u_steady, p_val)
    Jₚ = NonlinearSolveBase.nonlinearsolve_∂f_∂p(primal_prob, primal_prob.f, u_steady, p_val)

    # 3. Solve Jᵤ z = -Jₚ via the user-selected `alg.linsolve`. LinearSolve's
    #    `LinearProblem` machinery routes Float64 sparse through the chosen
    #    factorisation (default UMFPACK) and Dual sparse through
    #    `DualLinearCache` splitting (see `ift_solve` below).
    z = ift_solve(Jᵤ, -Jₚ, alg.linsolve)

    # 4. Chain rule: partials(u*[i]) = Σⱼ z[i,j] · partials(prob.p[j]).
    partials = [
        sum(z[i, j] * ForwardDiff.partials(prob.p[j]) for j in eachindex(prob.p))
            for i in eachindex(u_steady)
    ]

    # 5. Repackage u_value + partials as Dual{T,V,P} matching prob.p's flavour.
    dual_u = NonlinearSolveBase.nonlinearsolve_dual_solution(u_steady, partials, prob.p)

    return SciMLBase.build_solution(prob, alg, dual_u, sol.resid; retcode = sol.retcode)
end

# IFT solve dispatch.
#
# Both Float64 and Dual sparse `A` go through the same LinearSolve path so the
# user's `alg.linsolve` choice (set on `CTKAlg`) is honoured end-to-end —
# primal Newton iterations and the IFT back-solve use the same factorisation
# kind. Float64 sparse routes to the chosen factorisation directly; Dual
# sparse routes through LinearSolve's DualLinearCache, which splits
# A = A₀ + Σ Aₖεₖ, factorises A₀ once at value eltype, and back-solves for
# each partial direction reusing the factorisation — sparse all the way
# down, no densification.
#
# Looped per column of B because LinearSolve's `LinearProblem` takes a vector
# RHS; the factorisation in `cache` is reused across columns via `cache.b = …`.
function ift_solve(A, B, linsolve)
    Z = similar(B, promote_type(eltype(A), eltype(B)))
    cache = init(LinearProblem(A, copy(B[:, 1])), linsolve)
    Z[:, 1] = solve!(cache).u
    for j in 2:size(B, 2)
        cache.b = copy(B[:, j])
        Z[:, j] = solve!(cache).u
    end
    return Z
end

# ── NonlinearSolve / LinearSolve extension hooks ──────────────────────────────
# Bodies live in ext/AIBECSNonlinearSolveExt.jl and only become callable when
# the user has loaded NonlinearSolve and LinearSolve.
"""
    AIBECS.nonlinearproblem(prob::AbstractSteadyStateProblem) -> NonlinearProblem

Convert an AIBECS `SteadyStateProblem` into a `NonlinearProblem` with the
sparse Jacobian buffer type attached as `jac_prototype`. Without this, the
default Jacobian buffer is dense `Matrix{Float64}` and sparse linear-solver
factorisations (UMFPACK, KLU) error on the type mismatch. Requires
NonlinearSolve.jl and LinearSolve.jl to be loaded.
"""
function nonlinearproblem end

"""
    AIBECS.recommended_nlalg(; linsolve=UMFPACKFactorization())

Recommended NonlinearSolve algorithm for AIBECS-style sparse steady-state
problems: `NewtonRaphson` with `UMFPACKFactorization` (matches what
`LinearAlgebra.factorize` picks for sparse matrices today).
"""
function recommended_nlalg end
