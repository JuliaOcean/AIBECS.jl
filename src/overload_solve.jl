
struct CTKAlg <: SciMLBase.AbstractSteadyStateAlgorithm end

"""
    solve(prob::SciMLBase.AbstractSteadyStateProblem,
          alg::CTKAlg;
          nrm=norm,
          τstop=1e12*365*24*60*60,
          preprint="",
          maxItNewton=50)

Solves `prob` using the modified C.T.Kelley Shamanskii algorithm.
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
    SciMLBase.build_solution(prob, alg, x_steady, resid; retcode=:Success)
end

export solve, SteadyStateProblem, CTKAlg
