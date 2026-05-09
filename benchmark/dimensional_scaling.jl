# Diagonal (D_x, D_f) scaling of an AIBECS steady-state / nonlinear problem.
#
# Reference: Dennis, J. E., & Schnabel, R. B. (1996), "Numerical Methods for
# Unconstrained Optimization and Nonlinear Equations", SIAM Classics in Applied
# Mathematics 16; DOI 10.1137/1.9781611971200, Chapter 7 ("Stopping, Scaling,
# and Testing"). Newton iterates are mathematically unchanged under non-
# singular diagonal scaling (affine invariance, Deuflhard 2004 DOI
# 10.1007/978-3-642-23899-4); the convergence test on `‖F̃‖ ≤ tol` becomes a
# per-component relative-drift bound on the original problem.
#
# After `nondimensionalize`, both `ũ` and `F̃(ũ)` are O(1) per component, so
# every solver's default termination test ("‖F̃‖ ≤ abstol") becomes meaningful
# on a multi-tracer problem regardless of how unit-mismatched the original
# tracers are.

using SparseArrays
using LinearAlgebra
using SciMLBase

struct ScaledProblem{P, Su, SF, PS}
    original::P              # original SteadyStateProblem / NonlinearProblem
    scale_u::Su              # AbstractVector — per-component magnitude of u
    scale_F::SF              # AbstractVector — per-component magnitude of F
    scaled::PS               # rescaled problem of same kind, ũ-space
end

"""
    nondimensionalize(prob; scale_u, scale_F) -> ScaledProblem

Return a `ScaledProblem` wrapping `prob` (a `SteadyStateProblem` or
`NonlinearProblem`) under the change of variables ``ũ = u ./ scale_u``,
``F̃(ũ) = F(ũ .* scale_u, p) ./ scale_F``. The wrapped Jacobian is the
corresponding row-and-column-scaled `SparseMatrixCSC` (sparsity pattern
preserved).
"""
function nondimensionalize(prob; scale_u, scale_F)
    @assert length(scale_u) == length(prob.u0)
    @assert length(scale_F) == length(prob.u0)

    f_orig   = prob.f.f
    jac_orig = prob.f.jac
    inv_sF   = Diagonal(inv.(scale_F))
    Du       = Diagonal(scale_u)

    f_scaled(ũ, p, t = 0) = f_orig(ũ .* scale_u, p) ./ scale_F

    # In-place left/right Diagonal scaling on the freshly-computed sparse
    # Jacobian. `lmul!`/`rmul!` of a `Diagonal` against a `SparseMatrixCSC`
    # are O(nnz) and allocation-free (Julia 1.10+ stdlib), and preserve the
    # `SparseMatrixCSC{Float64,Int64}` type that UMFPACK/KLU expect.
    function jac_scaled(ũ, p, t = 0)
        J = jac_orig(ũ .* scale_u, p)
        lmul!(inv_sF, J)
        rmul!(J, Du)
        return J
    end

    ũ₀     = prob.u0 ./ scale_u
    f_new  = SciMLBase.ODEFunction(f_scaled; jac = jac_scaled)
    # `SciMLBase.remake` preserves the problem kind: a SteadyStateProblem in,
    # a SteadyStateProblem out — keeping CTKAlg's `solve(::SteadyStateProblem,
    # ::CTKAlg)` dispatch happy. Subsequent `AIBECS.nonlinearproblem(scaled)`
    # converts to a NonlinearProblem for the NewtonRaphson / wrapper rows.
    scaled = SciMLBase.remake(prob; f = f_new, u0 = ũ₀)

    return ScaledProblem(prob, scale_u, scale_F, scaled)
end

"""
    unscale(sp::ScaledProblem, ũ::AbstractVector) -> u

Round-trip a solution from scaled space back to the original problem's units.
"""
unscale(sp::ScaledProblem, ũ::AbstractVector) = ũ .* sp.scale_u

"""
    default_scales(prob) -> (scale_u, scale_F)

KINSOL-style heuristic (Hindmarsh et al. 2005): take per-component magnitudes
from the initial guess and the initial residual. Sensitive to a bad initial
guess but unit-agnostic and always works as a fallback for tracer setups that
don't supply a physics-based recipe.
"""
function default_scales(prob)
    F0 = prob.f.f(prob.u0, prob.p)
    return (
        max.(abs.(prob.u0), eps(eltype(prob.u0))),
        max.(abs.(F0),      eps(eltype(F0))),
    )
end

#= Quick perf check (run locally once; not in CI):
   J = ssprob.f.jac(ssprob.u0, ssprob.p)
   sF = Diagonal(rand(size(J,1))); sU = Diagonal(rand(size(J,2)))
   @benchmark $sF * $J * $sU
   @benchmark sparse($sF) * $J * sparse($sU)
   @benchmark (Jc = copy($J); lmul!($sF, Jc); rmul!(Jc, $sU))
   Expectation: in-place wins by 10-100×, allocates nothing. If not, swap.
=#
