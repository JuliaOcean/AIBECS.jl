using AIBECS
using Test
using ForwardDiff
using FiniteDiff
using NonlinearSolve
using LinearSolve
using SciMLBase

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

    # NL route: production `AIBECS.nonlinearproblem` attaches a *sparse*
    # `jac_prototype`, which currently triggers a bug in
    # `LinearSolveForwardDiffExt` on the nested-Dual (Hessian) path
    # (`MethodError: map!(partial_vals, ::Nothing, ::Vector{Float64})`).
    # The bug only fires for sparse Float64 prototype + Dual params; with
    # *dense* or no `jac_prototype` SciML's IFT works through both gradient
    # and Hessian. For this toy 24-box problem dense is fine — densifying
    # exercises the IFT logic end-to-end without the upstream sparse bug.
    function obj_nl(p)
        ssprob = SteadyStateProblem(F, x0, p)
        nlprob = SciMLBase.NonlinearProblem(ssprob)
        f_dense = SciMLBase.NonlinearFunction(
            nlprob.f.f;
            jac = nlprob.f.jac,
            jac_prototype = Matrix{Float64}(undef, length(x0), length(x0)),
        )
        prob = SciMLBase.NonlinearProblem(f_dense, nlprob.u0, nlprob.p)
        return sum(solve(prob, NewtonRaphson(); abstol = 1e-12).u)
    end

    # FD reference (uses obj_ctk; obj_nl returns the same primal value).
    g_fd = FiniteDiff.finite_difference_gradient(obj_ctk, p_ref)
    H_fd = FiniteDiff.finite_difference_hessian(obj_ctk, p_ref)

    # --- CTKAlg path: AIBECS's IFT dispatch added in src/overload_solve.jl ---
    g_ctk = ForwardDiff.gradient(obj_ctk, p_ref)
    H_ctk = ForwardDiff.hessian(obj_ctk, p_ref)            # nested Duals → IFT recursion
    @test g_ctk ≈ g_fd rtol = 1e-5
    @test H_ctk ≈ H_fd rtol = 1e-3

    # --- NonlinearSolve path: SciML's IFT dispatch (dense jac_prototype) ---
    g_nl = ForwardDiff.gradient(obj_nl, p_ref)
    H_nl = ForwardDiff.hessian(obj_nl, p_ref)
    @test g_nl ≈ g_fd rtol = 1e-5
    @test H_nl ≈ H_fd rtol = 1e-3

    # --- Cross-check: both AD paths use exact IFT, so they should agree
    #     much more tightly than FD does. ---
    @test g_ctk ≈ g_nl rtol = 1e-10
    @test H_ctk ≈ H_nl rtol = 1e-8

    # --- Sparse NL gradient: the production path through
    #     `AIBECS.nonlinearproblem` (sparse jac_prototype) works for the
    #     single-Dual gradient case. The nested-Dual Hessian via this same
    #     path is currently broken upstream in LinearSolve (see comment on
    #     `obj_nl` above) and is therefore not asserted here. ---
    obj_nl_sparse(p) = sum(
        solve(AIBECS.nonlinearproblem(SteadyStateProblem(F, x0, p)),
            AIBECS.recommended_nlalg(); abstol = 1e-12).u
    )
    g_nl_sparse = ForwardDiff.gradient(obj_nl_sparse, p_ref)
    @test g_nl_sparse ≈ g_fd rtol = 1e-5
    @test g_nl_sparse ≈ g_ctk rtol = 1e-10
end
