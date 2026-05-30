module demo

using AIBECS

struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end

struct SteadyStateSolution
    grd::OceanGrid
    prob::SteadyStateProblem
    state::AIBECS.SciMLBase.NonlinearSolution
end

function demo1()
    grd, TOCIM2 = OCIM2.load()

    z = depthvec(grd)
    function G(x, p)
        @unpack τ, z₀ = p
        return @. 1 - x / τ * (z ≤ z₀)
    end

    nb = sum(iswet(grd))  # number of wet boxes
    x_init = zeros(nb)    # Start with age = 0 everywhere

    p = IdealAgeParameters(1.0, z[1])
    F = AIBECSFunction(TOCIM2, G)

    prob = SteadyStateProblem(F, x_init, p)
    age = solve(prob, CTKAlg())
    age_in_yrs = age * u"s" .|> u"yr"

    SteadyStateSolution(grd,prob,age)
end

end

