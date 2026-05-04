
# Algorithms tested on every toy circulation. CTKAlg is the AIBECS-native
# solver and works without optional deps. The NonlinearSolve algorithms
# require the AIBECSNonlinearSolveExt extension to be loaded.
ctk_alg() = CTKAlg()
nlsolve_default()           = AIBECS.recommended_nlalg()
nlsolve_newton_umfpack()    = NewtonRaphson(linsolve = UMFPACKFactorization())
nlsolve_newton_klu()        = NewtonRaphson(linsolve = KLUFactorization())
nlsolve_trustregion_umfpack() = TrustRegion(linsolve = UMFPACKFactorization())
nlsolve_lm_umfpack()        = LevenbergMarquardt(linsolve = UMFPACKFactorization())

solver_cases = [
    (label = "CTKAlg",                       build = ctk_alg,                  needs_nlprob = false),
    (label = "recommended_nlalg",            build = nlsolve_default,          needs_nlprob = true),
    (label = "NewtonRaphson + UMFPACK",      build = nlsolve_newton_umfpack,   needs_nlprob = true),
    (label = "NewtonRaphson + KLU",          build = nlsolve_newton_klu,       needs_nlprob = true),
    (label = "TrustRegion + UMFPACK",        build = nlsolve_trustregion_umfpack, needs_nlprob = true),
    (label = "LevenbergMarquardt + UMFPACK", build = nlsolve_lm_umfpack,       needs_nlprob = true),
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
