
# Algorithms tested on every toy circulation. CTKAlg is the AIBECS-native
# solver and works without optional deps. The NonlinearSolve algorithms
# require the AIBECSNonlinearSolveExt extension to be loaded.
solver_cases = [
    (label = "CTKAlg",                  build = () -> CTKAlg(),                                       needs_nlprob = false),
    (label = "recommended_nlalg",       build = AIBECS.recommended_nlalg,                             needs_nlprob = true),
    (label = "NewtonRaphson + UMFPACK", build = () -> NewtonRaphson(linsolve = UMFPACKFactorization()), needs_nlprob = true),
    (label = "NewtonRaphson + KLU",     build = () -> NewtonRaphson(linsolve = KLUFactorization()),    needs_nlprob = true),
]

@testset "Solvers" begin
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n)
    @testset "$(case.label)" for case in solver_cases
        @testset "$(nx)x, $(np)p" for nx in 1:2, np in 1:2
            testp = AIBECS.reconstruct(TestParameters, np * vec(p))
            ssprob = SteadyStateProblem(fun, nx * x, testp)
            prob = case.needs_nlprob ? AIBECS.nonlinearproblem(ssprob) : ssprob
            s = solve(prob, case.build())
            @test s isa SciMLBase.AbstractSciMLSolution
            @test norm(s.u) / norm(fun.f(s.u, testp, 0)) > ustrip(upreferred(1e5u"Myr"))
        end
    end
end
