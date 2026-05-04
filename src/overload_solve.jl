
"""
    CTKAlg <: SciMLBase.AbstractSteadyStateAlgorithm

Marker algorithm that selects AIBECS's customised Newton–Chord–Shamanskii
solver when passed to `SciMLBase.solve(::SteadyStateProblem, ::CTKAlg)`.

```julia
prob = SteadyStateProblem(fun, x₀, p)
sol  = solve(prob, CTKAlg(); maxItNewton=50)
```

The underlying engine is `AIBECS.NewtonChordShamanskii` (with Armijo line
search and lazy Jacobian refresh, after Kelley 2003).
"""
struct CTKAlg <: SciMLBase.AbstractSteadyStateAlgorithm end

"""
    solve(prob::SciMLBase.AbstractSteadyStateProblem,
          alg::CTKAlg;
          nrm=norm,
          τstop=1e12*365*24*60*60,
          preprint="",
          maxItNewton=50)

Solves `prob` using an AIBECS-customized version of C.T.Kelley Shamanskii-method algorithm.
"""
function SciMLBase.solve(prob::SciMLBase.AbstractSteadyStateProblem,
                         alg::CTKAlg;
                         nrm=norm,
                         τstop=ustrip(u"s", 1e6u"Myr"),
                         preprint="",
                         maxItNewton=50)
    # Define the functions according to SciMLBase.SteadyStateProblem type
    p = prob.p
    x0 = copy(prob.u0)
    dx = copy(x0)
    function F(x)
        if SciMLBase.isinplace(prob)
            prob.f(dx, x, p)
            dx
        else
            prob.f(x, p)
        end
    end
    ∇ₓF(x) = prob.f.jac(x, p)
    # Compute `u_steady` and `resid` as per SciMLBase using my algorithm
    x_steady = NewtonChordShamanskii(F, ∇ₓF, nrm, x0, τstop; preprint=preprint, maxItNewton=maxItNewton)
    resid = F(x_steady)
    # Return the common SciMLBase solution type
    SciMLBase.build_solution(prob, alg, x_steady, resid; retcode=SciMLBase.ReturnCode.Success)
end

export solve, SteadyStateProblem, CTKAlg

# ── NonlinearSolve / LinearSolve extension hooks ──────────────────────────────
# Bodies live in ext/AIBECSNonlinearSolveExt.jl and only become callable when
# the user has loaded NonlinearSolve, LinearSolve, ADTypes, DifferentiationInterface,
# SparseConnectivityTracer, and SparseMatrixColorings.
"""
    AIBECS.nonlinearproblem(prob::AbstractSteadyStateProblem) -> NonlinearProblem

Build a `NonlinearProblem` from an AIBECS `SteadyStateProblem` with the
Jacobian sparsity prototype attached so LinearSolve can dispatch to a
sparse-aware factorization. Requires NonlinearSolve.jl to be loaded.
"""
function nonlinearproblem end

"""
    AIBECS.default_sparse_ad()

Default sparse automatic-differentiation backend used as a fallback when
no analytical Jacobian is supplied to NonlinearSolve. Combines ForwardDiff
(via DifferentiationInterface), SparseConnectivityTracer for sparsity
detection, and SparseMatrixColorings for graph coloring.
"""
function default_sparse_ad end

"""
    AIBECS.recommended_nlalg(; linsolve=UMFPACKFactorization(), autodiff=default_sparse_ad())

Recommended NonlinearSolve algorithm for AIBECS-style sparse steady-state
problems: `NewtonRaphson` with `UMFPACKFactorization` (matches what
`LinearAlgebra.factorize` picks for sparse matrices today). `autodiff` is
only consulted if the supplied `prob.f.jac` is dropped.
"""
function recommended_nlalg end
