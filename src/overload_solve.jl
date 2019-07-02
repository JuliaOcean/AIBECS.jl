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
                          τstop=ustrip(1e6u"Myr"|>u"s"),
                          preprint="",
                          maxItNewton=50)
    # Define the functions according to DiffEqBase.SteadyStateProblem type
    p = prob.p
    t = 0
    x0 = copy(prob.u0)
    dx, df = copy(x0), copy(x0)
    F(x) = prob.f(dx, x, p, t)
    ∇ₓF(x) = prob.f(df, dx, x, p, t)
    # Compute `u_steady` and `resid` as per DiffEqBase using my algorithm
    x_steady = NewtonChordShamanskii(F, ∇ₓF, nrm, x0, τstop; preprint=preprint, maxItNewton=maxItNewton)
    resid = F(x_steady)
    # Return the common DiffEqBase solution type
    DiffEqBase.build_solution(prob, alg, x_steady, resid; retcode=:Success)
end

"""
    SteadyStateProblem(F, ∇ₓF, x, p)

Returns the `SteadyStateProblem` defined by `F(x,p)=0`.
"""
function DiffEqBase.SteadyStateProblem(F, ∇ₓF, x, p)
    f(dx, x, p, t) = F(x, p)
    f(df, dx, x, p, t) = ∇ₓF(x, p)
    return DiffEqBase.SteadyStateProblem(f, x, p)
end

export solve, SteadyStateProblem, CTKAlg, CTKAlg2
