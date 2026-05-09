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
CTKAlg(; linsolve = nothing) = CTKAlg(linsolve)

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
    # Default linsolve mirrors `LinearAlgebra.factorize` on a sparse Float64
    # Jacobian (which lands on UMFPACK).
    linsolve_alg = alg.linsolve === nothing ? UMFPACKFactorization() : alg.linsolve
    # Compute `u_steady` and `resid` as per SciMLBase using my algorithm
    x_steady = NewtonChordShamanskii(
        F, ∇ₓF, x0, linsolve_alg, tc_cache;
        preprint = preprint, maxItNewton = maxItNewton
    )
    resid = F(x_steady)
    # Return the common SciMLBase solution type
    return SciMLBase.build_solution(prob, alg, x_steady, resid; retcode = SciMLBase.ReturnCode.Success)
end

export solve, SteadyStateProblem, CTKAlg

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
