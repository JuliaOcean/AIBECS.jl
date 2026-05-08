@info "Running test/solvers.jl — steady-state solvers" Circulation
# Algorithms tested on every toy circulation. CTKAlg is the AIBECS-native
# solver and works without optional deps. The NonlinearSolve algorithms
# require the AIBECSNonlinearSolveExt extension to be loaded. CTKAlg now
# accepts a `linsolve` kwarg routed through LinearSolve, so we exercise
# both the default (UMFPACK) and KLU linear-solver paths.
solver_cases = [
    (label = "CTKAlg",                  build = () -> CTKAlg(),                                 needs_nlprob = false),
    (label = "CTKAlg + KLU",            build = () -> CTKAlg(linsolve = KLUFactorization()),    needs_nlprob = false),
    (label = "recommended_nlalg",       build = AIBECS.recommended_nlalg,                       needs_nlprob = true),
    (label = "NewtonRaphson + UMFPACK", build = () -> NewtonRaphson(linsolve = UMFPACKFactorization()), needs_nlprob = true),
    (label = "NewtonRaphson + KLU",     build = () -> NewtonRaphson(linsolve = KLUFactorization()),     needs_nlprob = true),
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
            # Residual is small in absolute ∞-norm. NonlinearSolveBase's
            # default abstol for Float64 is ~eps^0.8 ≈ 7e-13; 1e-6 leaves
            # plenty of margin.
            @test maximum(abs, fun.f(s.u, testp, 0)) ≤ 1.0e-6
        end
    end

    # Smoke test: CTKAlg accepts a NonlinearSolveBase termination mode and
    # honors `reltol`. NormTerminationMode passes if either abstol or reltol
    # is satisfied; we set abstol = 0 to force the reltol path.
    @testset "CTKAlg termination_condition (reltol path)" begin
        testp = p
        ssprob = SteadyStateProblem(fun, x, testp)
        s = solve(
            ssprob, CTKAlg();
            termination_condition = NonlinearSolveBase.NormTerminationMode(norm),
            abstol = 0, reltol = 1.0e-8,
        )
        @test s isa SciMLBase.AbstractSciMLSolution
        residual = fun.f(s.u, testp, 0)
        @test norm(residual) ≤ 1.0e-8 * norm([residual; s.u])
    end
end
