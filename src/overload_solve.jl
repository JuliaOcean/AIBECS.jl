using DiffEqBase

struct CTKAlg <: DiffEqBase.AbstractSteadyStateAlgorithm end

"""
    solve(prob::DiffEqBase.AbstractSteadyStateProblem,
          alg::CTKAlg;
          nrm=norm,
          τstop=1e12*365*24*60*60,
          preprint="",
          maxItNewton=50)

Solves `prob` using the modified C.T.Kelley Shamanskii algorithm.
"""
function DiffEqBase.solve(prob::DiffEqBase.AbstractSteadyStateProblem,
                          alg::CTKAlg;
                          nrm=norm,
                          τstop=ustrip(u"s", 1e6u"Myr"),
                          preprint="",
                          maxItNewton=50)
    # Define the functions according to DiffEqBase.SteadyStateProblem type
    p = prob.p
    x0 = copy(prob.u0)
    dx = copy(x0)
    F(x) = prob.f.f(x, p, 0)
    ∇ₓF(x) = prob.f.jac(x, p, 0)
    # Compute `u_steady` and `resid` as per DiffEqBase using my algorithm
    x_steady = NewtonChordShamanskii(F, ∇ₓF, nrm, x0, τstop; preprint=preprint, maxItNewton=maxItNewton)
    resid = F(x_steady)
    # Return the common DiffEqBase solution type
    DiffEqBase.build_solution(prob, alg, x_steady, resid; retcode=:Success)
end

# Hacky way for F1Method to construct CTK-solvable SteadyStateProblems from F and ∇ₓF
# TODO add test for it (maybe a F1Method test in AIBECS? Circular but oh well...
function DiffEqBase.SteadyStateProblem(F, ∇ₓF, x, p)
    f(u, p, t) = F(u, p)
    jac(u, p, t) = ∇ₓF(u, p)
    SteadyStateProblem(ODEFunction(f, jac=jac), x, p)
end

export solve, SteadyStateProblem, CTKAlg
