using ForwardDiff: Dual

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

"""
    solve(prob::SciMLBase.AbstractSteadyStateProblem,
          alg::CTKAlg;
          termination_condition = AbsNormSafeBestTerminationMode(Base.Fix1(maximum, abs)),
          abstol = nothing,
          reltol = nothing,
          preprint = "",
          maxItNewton = 50)

Solves `prob` using an AIBECS-customized Newton–Chord–Shamanskii method with
Armijo line search. Convergence is delegated to NonlinearSolveBase termination
modes — the default mirrors `NonlinearSolve`'s `:regular` algorithm default
(`AbsNormSafeBestTerminationMode` on the ∞-norm of the residual). Pass
`abstol`/`reltol` to tighten or loosen, or pass any
`AbstractNonlinearTerminationMode` (e.g. `NormTerminationMode(norm)`,
`RelNormTerminationMode(norm)`) to swap criteria.
"""
function SciMLBase.solve(
        prob::SciMLBase.AbstractSteadyStateProblem,
        alg::CTKAlg;
        termination_condition = NonlinearSolveBase.AbsNormSafeBestTerminationMode(
            Base.Fix1(maximum, abs)
        ),
        abstol = nothing,
        reltol = nothing,
        preprint = "",
        maxItNewton = 50
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
    x_steady = NewtonChordShamanskii(
        F, ∇ₓF, x0, alg.linsolve, tc_cache;
        preprint = preprint, maxItNewton = maxItNewton
    )
    resid = F(x_steady)
    # Return the common SciMLBase solution type
    return SciMLBase.build_solution(prob, alg, x_steady, resid; retcode = SciMLBase.ReturnCode.Success)
end

export solve, SteadyStateProblem, CTKAlg

# ── AD-through-the-solve via the implicit function theorem ───────────────────
# Mirror of `NonlinearSolveBase.nonlinearsolve_forwarddiff_solve` (defined in
# `NonlinearSolveBaseForwardDiffExt`) for AIBECS's `CTKAlg` route. SciML's IFT
# dispatch only fires for `NonlinearProblem`/`NonlinearLeastSquaresProblem`;
# `SteadyStateProblem` with a Dual-typed `p` does not hit it, so we add the
# equivalent dispatch here.
const DualSteadyStateProblem = SciMLBase.SteadyStateProblem{
    <:Union{Number, <:AbstractArray}, iip,
    <:Union{<:Dual{T, V, P}, <:AbstractArray{<:Dual{T, V, P}}},
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
"""
function SciMLBase.solve(prob::DualSteadyStateProblem, alg::CTKAlg, args...; kwargs...)
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

    # 3. Solve Jᵤ z = -Jₚ. Multiple-dispatch on eltype routes Float64 sparse
    #    through UMFPACK (one-shot, not a hot path — fine to bypass the
    #    `alg.linsolve` indirection) and Dual eltype through dense generic LU
    #    (UMFPACK is Float64-only, and SparseMatrixCSC{Dual} has no sparse LU).
    z = ift_solve(Jᵤ, -Jₚ)

    # 4. Chain rule: partials(u*[i]) = Σⱼ z[i,j] · partials(prob.p[j]).
    partials = [
        sum(z[i, j] * ForwardDiff.partials(prob.p[j]) for j in eachindex(prob.p))
            for i in eachindex(u_steady)
    ]

    # 5. Repackage u_value + partials as Dual{T,V,P} matching prob.p's flavour.
    dual_u = NonlinearSolveBase.nonlinearsolve_dual_solution(u_steady, partials, prob.p)

    return SciMLBase.build_solution(prob, alg, dual_u, sol.resid; retcode = sol.retcode)
end

# IFT solve dispatch. Float64 sparse → UMFPACK via `\`. Dual sparse →
# densify and use Julia's generic dense LU (UMFPACK is Float64-only and
# `SparseMatrixCSC{<:Dual}` has no sparse LU implementation, so `\` would
# infinite-recurse). Densification is O(n²) memory and only hit on the nested
# -Dual Hessian path, which already pays an `O(n)`-partials overhead — this is
# acceptable for the small/medium circulations where Hessian-via-AD is even
# tractable; production-scale Hessians should use F1Method's analytical formula.
ift_solve(A, B) = A \ B
ift_solve(A::AbstractSparseMatrix{<:Dual}, B) = Matrix(A) \ B

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
