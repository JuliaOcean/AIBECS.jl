module AIBECSNonlinearSolveExt

using AIBECS
using SciMLBase
using NonlinearSolve
using LinearSolve
using ADTypes
using DifferentiationInterface
using SparseConnectivityTracer
using SparseMatrixColorings
using ForwardDiff

function AIBECS.nonlinearproblem(prob::SciMLBase.AbstractSteadyStateProblem)
    # TODO(NonlinearSolve interop): AIBECS's f and jac are 3-arg ODE-style
    # (u, p, t=0). NonlinearSolve's AutoSpecialize wraps the function via
    # FunctionWrappers and re-promotes to in-place internally
    # (NonlinearSolveBase/src/autospecialize.jl:wrapfun_iip), regardless of
    # the SteadyStateProblem's iip flag. Calling our 3-arg f through that
    # in-place wrapper shifts (u,p,t) → (du,u,p) and the parameters object
    # ends up in the wrong slot, blowing up inside @unpack. Explicitly
    # wrapping in 2-arg shims and pinning NonlinearFunction{false} /
    # NonlinearProblem{false} sidesteps the iip-promotion path.
    f_ode   = prob.f.f
    jac_ode = prob.f.jac
    f_nl(u, p)   = f_ode(u, p)
    jac_nl(u, p) = jac_ode(u, p)
    proto = jac_ode(prob.u0, prob.p)
    nf = SciMLBase.NonlinearFunction{false}(f_nl; jac = jac_nl, jac_prototype = proto)
    return SciMLBase.NonlinearProblem{false}(nf, prob.u0, prob.p)
end

AIBECS.default_sparse_ad() = AutoSparse(
    AutoForwardDiff();
    sparsity_detector = TracerSparsityDetector(),
    coloring_algorithm = GreedyColoringAlgorithm(),
)

AIBECS.recommended_nlalg(;
    linsolve = UMFPACKFactorization(),
    autodiff = AIBECS.default_sparse_ad(),
) = NewtonRaphson(; linsolve, autodiff)

end # module
