module AIBECSNonlinearSolveExt

using AIBECS
using SciMLBase
using NonlinearSolve
using LinearSolve

# Convert an AIBECS SteadyStateProblem to a NonlinearProblem and attach the
# Jacobian sparsity *type* via `jac_prototype`. The latter is what tells
# NonlinearSolve's linear-solve cache to allocate a SparseMatrixCSC buffer for
# J — without it the buffer defaults to dense Matrix{Float64}, and sparse
# factorisations like UMFPACK/KLU then error trying to read `.nzval` off it.
function AIBECS.nonlinearproblem(prob::SciMLBase.AbstractSteadyStateProblem)
    nlprob = SciMLBase.NonlinearProblem(prob)
    jac_prototype = prob.f.jac(prob.u0, prob.p)
    f = SciMLBase.NonlinearFunction(nlprob.f.f; jac = nlprob.f.jac, jac_prototype)
    return SciMLBase.NonlinearProblem(f, nlprob.u0, nlprob.p)
end

AIBECS.recommended_nlalg(; linsolve = UMFPACKFactorization()) = NewtonRaphson(; linsolve)

end # module
