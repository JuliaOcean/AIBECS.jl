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
    t = 0
    x0 = copy(prob.u0)
    dx = copy(x0)
    F(x) = prob.f.f(dx, x, p, t)
    ∇ₓF(x) = prob.f.jac(x, p, t)
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
function DiffEqBase.SteadyStateProblem(F, ∇ₓF, x, p::AbstractParameters)
    return SteadyStateProblem(ODEFunction(F, ∇ₓF, x, p), x, p)
end

function DiffEqBase.ODEProblem(F, ∇ₓF, x, p::AbstractParameters; tspan=(0.0, ustrip(u"s", 1000.0u"yr")))
    return ODEProblem(ODEFunction(F, ∇ₓF, x, p), x, tspan, p)
end

function DiffEqBase.ODEFunction(F, ∇ₓF, x, p::AbstractParameters)
    f = (du,u,p,t) -> du .= F(u,p)
    jac = (u,p,t) -> ∇ₓF(u,p)
    jp = ∇ₓF(x,p)
    return ODEFunction(f; jac=jac, jac_prototype=jp)
end

export solve, SteadyStateProblem, CTKAlg
