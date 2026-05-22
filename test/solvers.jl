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

    @testset "CTKAlg retcode" begin
        testp = p
        ssprob = SteadyStateProblem(fun, x, testp)

        s_ok = solve(ssprob, CTKAlg())
        @test SciMLBase.successful_retcode(s_ok.retcode)

        s_short = solve(ssprob, CTKAlg(); maxItNewton = 1)
        @test s_short.retcode == SciMLBase.ReturnCode.MaxIters
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

@testset "Time-stepping" begin
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x0 = xgeo * ones(n)
    testp = p

    # Reference steady state via CTKAlg.
    ssprob = SteadyStateProblem(fun, x0, testp)
    u_ss = solve(ssprob, CTKAlg()).u

    # Long-time integration window: ~10× the slowest model timescale (τgeo
    # = 1 Myr) is enough that ImplicitEuler's L-stable damping reaches the
    # steady state to ~percent relative accuracy on toy circulations.
    t_end = ustrip(upreferred(10.0u"Myr"))
    dt    = ustrip(upreferred(100.0u"kyr"))

    @testset "ODEProblem → AIBECS.odeproblem attaches sparse jac_prototype" begin
        odeprob = AIBECS.odeproblem(fun, x0, (0.0, t_end), testp)
        @test odeprob isa SciMLBase.ODEProblem
        @test odeprob.f.jac_prototype isa SparseMatrixCSC
    end

    @testset "ImplicitEuler long-time → CTKAlg steady state" begin
        odeprob = AIBECS.odeproblem(fun, x0, (0.0, t_end), testp)
        sol = solve(
            odeprob, ImplicitEuler(linsolve = UMFPACKFactorization());
            dt = dt, saveat = [t_end],
        )
        @test SciMLBase.successful_retcode(sol.retcode)
        @test maximum(abs, sol.u[end] - u_ss) ≤ 1.0e-2 * maximum(abs, u_ss)
    end

    @testset "splitodeproblem jac is constant in u" begin
        # Use the bgc L/NL split to confirm the structural promise:
        # the linear part's Jacobian (-T + ∇ₓL) does not depend on u.
        splitfun = AIBECSSplitFunction(T_all, L_all, NL_all, nb)
        sprob = AIBECS.splitodeproblem(splitfun, x0, (0.0, t_end), testp)
        x_alt = x0 .* (1.0 .+ 0.1 .* randn(n))
        @test sprob.f.jac(x0, testp, 0.0) === sprob.f.jac(x_alt, testp, 0.0)
    end

    # IMEX cross-check on a linear-only model.
    #
    # The bgc test model puts Michaelis-Menten uptake (stiff, τDIP ≈ 30 days)
    # into the *nonlinear* part of the split — which IMEX schemes like SBDF2
    # treat explicitly. On a 100-kyr step that's catastrophically unstable
    # (the explicit dt would have to be a fraction of τDIP). To verify the
    # IMEX wiring honestly, we build an inline single-tracer linear model
    # where all reactive dynamics live in L and NL is zero. SBDF2 then
    # reduces to a stable implicit-multistep march, and we can cross-check
    # against the same circulation's steady-state solver.
    @testset "SBDF2 (IMEX) on linear model → CTKAlg steady state" begin
        τ_lin = ustrip(upreferred(1.0u"yr"))
        Llin(x, p) = (xgeo .- x) / τ_lin
        NLzero(x, p) = zero(x)
        Glin(x, p) = Llin(x, p)
        funL = AIBECSFunction(T_D, Glin, nb)
        splitfunL = AIBECSSplitFunction(T_D, Llin, NLzero, nb)
        x0L = xgeo * ones(nb)
        u_ss_L = solve(SteadyStateProblem(funL, x0L, testp), CTKAlg()).u
        sprobL = AIBECS.splitodeproblem(splitfunL, x0L, (0.0, t_end), testp)
        sol = solve(
            sprobL, SBDF2(linsolve = UMFPACKFactorization());
            dt = dt, saveat = [t_end],
        )
        @test SciMLBase.successful_retcode(sol.retcode)
        @test maximum(abs, sol.u[end] - u_ss_L) ≤ 1.0e-2 * maximum(abs, u_ss_L)
    end
end
