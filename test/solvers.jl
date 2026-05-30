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
            # Residual is small in absolute ∞-norm. CTKAlg's default
            # `NormTerminationMode` with `reltol ≈ 1/(1 Myr in seconds)`
            # ≈ 3.17e-14 converges well below 1e-6 for these toy
            # problems; NonlinearSolve's `:regular` default (~eps^0.8
            # ≈ 7e-13 absolute) is similarly loose against this margin.
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

    @testset "CTKAlg retcode" begin
        testp = p
        ssprob = SteadyStateProblem(fun, x, testp)

        # Success path: no @warn, successful retcode.
        s_ok = @test_logs solve(ssprob, CTKAlg())
        @test SciMLBase.successful_retcode(s_ok.retcode)

        # MaxIters path: @warn fires, retcode = MaxIters.
        s_short = @test_logs (:warn, r"CTKAlg did not converge") solve(
            ssprob, CTKAlg(); maxItNewton = 1,
        )
        @test s_short.retcode == SciMLBase.ReturnCode.MaxIters
    end

    @testset "CTKAlg stats" begin
        testp = p
        ssprob = SteadyStateProblem(fun, x, testp)
        s = solve(ssprob, CTKAlg())
        @test s.stats isa SciMLBase.NLStats
        @test s.stats.nsteps > 0
        # Documented mapping: every Jacobian refresh refactors with the
        # direct-LU linsolves used here, so `nfactors == njacs`.
        @test s.stats.nfactors == s.stats.njacs
    end

    @testset "default_termination_condition" begin
        dtc = default_termination_condition()
        @test keys(dtc) == (:termination_condition, :abstol, :reltol)
        @test dtc.abstol == 0.0
        @test dtc.reltol > 0
        # Splattable into solve(...) without erroring.
        ssprob = SteadyStateProblem(fun, x, p)
        s = solve(ssprob, CTKAlg(); dtc...)
        @test SciMLBase.successful_retcode(s.retcode)
    end

    @testset "CTKAlg warm-start via linear_cache" begin
        testp = p
        ssprob = SteadyStateProblem(fun, x, testp)

        # Cold solve as reference.
        s_cold = solve(ssprob, CTKAlg())

        # Build a vector-RHS LinearCache at the converged Jacobian and force
        # factorisation. After solve!, isfresh = false and cacheval holds the
        # UMFPACK factors — exactly the state a warm-start caller would set up.
        lc = init(
            LinearProblem(fun.jac(s_cold.u, testp), similar(s_cold.u)),
            UMFPACKFactorization(),
        )
        solve!(lc)

        # Same-point warm-start: residual is similarly small. State agreement
        # is bounded by the solver tolerance (~1e-6 here) — both Newton runs
        # terminate at a finite residual norm, not at machine precision, so
        # iterates can drift at that scale.
        s_warm = solve(ssprob, CTKAlg(); linear_cache = lc)
        @test maximum(abs, fun.f(s_warm.u, testp, 0)) ≤ 1.0e-6
        @test maximum(abs, s_warm.u - s_cold.u) ≤ 1.0e-5

        # Nearby-parameter warm-start: re-solve at a perturbed p with the
        # same (now slightly stale) factors. Shamanskii will refresh if the
        # ratio check fires; either way convergence must still happen.
        testp2 = AIBECS.reconstruct(TestParameters, 1.05 * vec(p))
        ssprob2 = SteadyStateProblem(fun, s_cold.u, testp2)
        s_warm2 = solve(ssprob2, CTKAlg(); linear_cache = lc)
        @test maximum(abs, fun.f(s_warm2.u, testp2, 0)) ≤ 1.0e-6
    end
end
