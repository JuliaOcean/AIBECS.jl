module AIBECSNonlinearSolveExt

using AIBECS
using SciMLBase
using NonlinearSolve
using NonlinearSolveBase
using LinearSolve

# Convert an AIBECS SteadyStateProblem to a NonlinearProblem and attach the
# Jacobian sparsity *type* via `jac_prototype`. The latter is what tells
# NonlinearSolve's linear-solve cache to allocate a SparseMatrixCSC buffer for
# J — without it the buffer defaults to dense Matrix{Float64}, and sparse
# factorisations like UMFPACK/KLU then error trying to read `.nzval` off it.
#
# `prob.p` may carry `ForwardDiff.Dual` numbers when this conversion is invoked
# inside `ForwardDiff.gradient`/`.hessian`. SciML's
# `__solve(::DualNonlinearProblem, …)` strips the duals and runs the *primal*
# solve at Float64, but the Jacobian buffer attached here is what flows into
# LinearSolve — UMFPACK/KLU reject Dual eltypes. Build the prototype from
# value-typed inputs so the primal path stays Float64. `nodual_value` is the
# identity for non-Dual inputs, so the existing non-AD code path is unaffected.
function AIBECS.nonlinearproblem(prob::SciMLBase.AbstractSteadyStateProblem)
    nlprob = SciMLBase.NonlinearProblem(prob)
    p_val = NonlinearSolveBase.nodual_value(prob.p)
    u0_val = NonlinearSolveBase.nodual_value(prob.u0)
    jac_prototype = prob.f.jac(u0_val, p_val)
    f = SciMLBase.NonlinearFunction(nlprob.f.f; jac = nlprob.f.jac, jac_prototype)
    return SciMLBase.NonlinearProblem(f, nlprob.u0, nlprob.p)
end

AIBECS.recommended_nlalg(; linsolve = UMFPACKFactorization()) = NewtonRaphson(; linsolve)

end # module
