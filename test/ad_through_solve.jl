using AIBECS
using Test
using ForwardDiff
using FiniteDiff
using NonlinearSolve
using LinearSolve

@testset "AD through solve (IFT)" begin
    grd, T = Primeau_2x2x2.load()
    nb = sum(iswet(grd))

    # Single-tracer damped restoring source. Damping at rate 1/p[1] is
    # applied to every wet box (not just the surface), so the steady state
    # is well-conditioned for `T + I/p[1]` regardless of T's null space.
    # Analytical ∂G keeps the inner Jacobian assembly free of nested Duals.
    G(x, p) = @. p[2] - x / p[1]
    ∂G(x, p) = @. -1 / p[1] * one(x)
    F = AIBECSFunction(T, G; ∂G = ∂G)
    x0 = zeros(nb)
    p_ref = [1.0, 1.0]

    obj_ctk(p) = sum(
        solve(SteadyStateProblem(F, x0, p), CTKAlg();
            abstol = 1e-12, maxItNewton = 200).u
    )

    # NL route — bare `NonlinearProblem(SteadyStateProblem(...))`, no
    # `jac_prototype` attached. `AIBECS.nonlinearproblem` would attach a *sparse*
    # `jac_prototype`, which still trips `LinearSolveForwardDiffExt.setu!` on
    # the nested-Dual (Hessian) path even after LinearSolve v3.76.0 fixed the
    # related `setb!` bug (SciML/LinearSolve.jl#972). For this 24-box problem
    # the no-prototype path costs nothing meaningful; production code that
    # wants the sparse prototype can still use `AIBECS.nonlinearproblem` for
    # the primal/gradient case (asserted separately below).
    obj_nl(p) = sum(
        solve(NonlinearProblem(SteadyStateProblem(F, x0, p)),
            NewtonRaphson(); abstol = 1e-12).u
    )

    # FD reference (uses obj_ctk; obj_nl returns the same primal value).
    g_fd = FiniteDiff.finite_difference_gradient(obj_ctk, p_ref)
    H_fd = FiniteDiff.finite_difference_hessian(obj_ctk, p_ref)

    # --- CTKAlg path: AIBECS's IFT dispatch in src/overload_solve.jl ---
    g_ctk = ForwardDiff.gradient(obj_ctk, p_ref)
    H_ctk = ForwardDiff.hessian(obj_ctk, p_ref)            # nested Duals → IFT recursion
    @test g_ctk ≈ g_fd rtol = 1e-5
    @test H_ctk ≈ H_fd rtol = 1e-3

    # --- NonlinearSolve path: SciML's IFT dispatch (sparse jac_prototype) ---
    g_nl = ForwardDiff.gradient(obj_nl, p_ref)
    H_nl = ForwardDiff.hessian(obj_nl, p_ref)
    @test g_nl ≈ g_fd rtol = 1e-5
    @test H_nl ≈ H_fd rtol = 1e-3

    # --- Cross-check: both AD paths use exact IFT, so they should agree
    #     much more tightly than FD does. ---
    @test g_ctk ≈ g_nl rtol = 1e-10
    @test H_ctk ≈ H_nl rtol = 1e-8

    # --- Sparse `AIBECS.nonlinearproblem` path: works for primal + gradient
    #     (single-Dual layer). Hessian via this same path still trips
    #     `LinearSolveForwardDiffExt.setu!` upstream and is not asserted here. ---
    obj_nl_sparse(p) = sum(
        solve(AIBECS.nonlinearproblem(SteadyStateProblem(F, x0, p)),
            AIBECS.recommended_nlalg(); abstol = 1e-12).u
    )
    g_nl_sparse = ForwardDiff.gradient(obj_nl_sparse, p_ref)
    @test g_nl_sparse ≈ g_fd rtol = 1e-5
    @test g_nl_sparse ≈ g_ctk rtol = 1e-10

    # --- IFT dispatch without analytical `∂Gs`: same setup but the inner
    #     Jacobian falls back to ForwardDiff. The CTKAlg IFT dispatch only
    #     touches `prob.f.jac`, so it should be agnostic to how that jac was
    #     built. Cross-checked against `g_ctk` (analytical) at tight rtol. ---
    F_ad = AIBECSFunction(T, G)  # no ∂G kwarg → ForwardDiff jac
    obj_ctk_ad(p) = sum(
        solve(SteadyStateProblem(F_ad, x0, p), CTKAlg();
            abstol = 1e-12, maxItNewton = 200).u
    )
    g_ctk_ad = ForwardDiff.gradient(obj_ctk_ad, p_ref)
    @test g_ctk_ad ≈ g_ctk rtol = 1e-10
end
