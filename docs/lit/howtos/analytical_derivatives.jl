#---------------------------------------------------------
# # [Supplying analytical derivatives for `G`](@id analytical_derivatives)
#---------------------------------------------------------

# By default, `AIBECSFunction` builds the Jacobian of the source–sink term
# `G` using ForwardDiff. For models where you already know the derivatives
# analytically, you can skip AD entirely by passing an `∂Gs` keyword. The
# analytical Jacobian is faster (no AD overhead per evaluation) and gives a
# cleaner path when you want to differentiate through `solve` itself, since
# it removes one level of AD nesting.

# ## API summary
#
# - Single-tracer call: pass `∂G = (x, p) -> diag_vec` to `AIBECSFunction`,
#   where `diag_vec` is the length-`nb` diagonal of `∂G/∂x`.
# - Multi-tracer call: pass `∂Gs` as a tuple-of-tuples mirroring `Gs`, where
#   `∂Gs[i][j](x₁, …, xₙₜ, p)` returns the diagonal of `∂Gᵢ/∂xⱼ`. **All
#   `nt × nt` blocks must be supplied** — use `zero(xⱼ)` for blocks that are
#   identically zero.
# - The diagonal-only assumption (each block is pointwise in space) is the
#   same as in the AD path, so structurally the resulting sparse Jacobian is
#   identical.
#
# Notation: `∇ₓG` is the assembled *Jacobian matrix*; `∂G` (and each
# `∂Gs[i][j]`) is a *diagonal vector* — the same partial derivative, packed
# differently.

# ## Example: ideal age, single tracer

using AIBECS
grd, T = Primeau_2x2x2.load()
z = depthvec(grd)

function G(x, p)
    @unpack τ, z₀ = p
    return @. 1 - x / τ * (z ≤ z₀)
end

# The derivative is pointwise, with `∂G/∂x = -1/τ` in the surface layer and
# `0` everywhere else:

function ∂G(x, p)
    @unpack τ, z₀ = p
    return @. -1 / τ * (z ≤ z₀) * one(x)
end

struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end

p = IdealAgeParameters(1.0, z[1])

# Construct the function with `∂G`, then solve as usual:

F = AIBECSFunction(T, G; ∂G = ∂G)
nb = sum(iswet(grd))
x_init = zeros(nb)
prob = SteadyStateProblem(F, x_init, p)
sol = solve(prob, CTKAlg())

# ## Validating analytical derivatives
#
# Hand-written derivatives are easy to get wrong. Use [`check_∂Gs`](@ref) to
# compare against ForwardDiff at a sample `(x, p)`. The returned `maxabs`
# should be at machine precision; the helper also returns both Jacobians as
# `ad_jac` and `analytical_jac` so you can inspect any disagreement
# block-by-block.

check_∂Gs(G, ∂G, sol.u, p, nb).maxabs

# ## Multi-tracer pattern
#
# For an `nt`-tracer model, `Gs = (G₁, G₂, …)` and `∂Gs` mirrors that shape:
#
# ```julia
# ∂Gs = (
#     (∂G1_∂x1, ∂G1_∂x2, …, ∂G1_∂xₙₜ),  # row 1: derivatives of G₁
#     (∂G2_∂x1, ∂G2_∂x2, …, ∂G2_∂xₙₜ),  # row 2: derivatives of G₂
#     …
# )
# fun = AIBECSFunction(Ts, Gs, nb; ∂Gs = ∂Gs)
# ```
#
# Each `∂Gᵢ_∂xⱼ` has the same call signature as `Gᵢ` — i.e.
# `(x₁, x₂, …, xₙₜ, p)` — even if the derivative does not depend on every
# tracer. Returning `zero(xⱼ)` is the idiom for blocks that are identically
# zero; this keeps shape uniform and allows derivative-of-derivative
# computations to flow through cleanly.
