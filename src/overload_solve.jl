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
    function F(x)
        if DiffEqBase.isinplace(prob)
            prob.f(dx, x, p)
            dx
        else
            prob.f(x, p)
        end
    end
    ∇ₓF(x) = prob.f(Val{:jac}, x, p)
    # Compute `u_steady` and `resid` as per DiffEqBase using my algorithm
    x_steady = NewtonChordShamanskii(F, ∇ₓF, nrm, x0, τstop; preprint=preprint, maxItNewton=maxItNewton)
    resid = F(x_steady)
    # Return the common DiffEqBase solution type
    DiffEqBase.build_solution(prob, alg, x_steady, resid; retcode=:Success)
end

export solve, SteadyStateProblem, CTKAlg
