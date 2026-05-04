#---------------------------------------------------------
# # [NonlinearSolve.jl integration](@id nonlinearsolve)
#---------------------------------------------------------

# AIBECS ships with its own steady-state solver, `CTKAlg()`, which has zero
# extra dependencies and is the default for every tutorial.
# For users who want to tap into the broader SciML ecosystem, AIBECS also
# provides an opt-in integration with [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl)
# and [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl).
# Both are *weak* dependencies: the integration only activates once you
# load them yourself.

using AIBECS
using NonlinearSolve, LinearSolve
using JLD2 # required by `OCIM0.load`

# We use the small OCIM0 circulation here so the example stays light;
# any AIBECS circulation works the same way.

grd, T_OCIM0 = OCIM0.load()

# Build the standard ideal-age problem (see the [Ideal age tutorial](@ref idealage)
# for the full derivation).

z = depthvec(grd)
function G(x, p)
    @unpack τ, z₀ = p
    return @. 1 - x / τ * (z ≤ z₀)
end

struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end

p = IdealAgeParameters(1.0, 30.0)
F = AIBECSFunction(T_OCIM0, G)
nb = sum(iswet(grd))
x_init = zeros(nb)
prob = SteadyStateProblem(F, x_init, p)

# ## Solving with NonlinearSolve

# Wrap the `SteadyStateProblem` with `AIBECS.nonlinearproblem` before
# calling `solve`. The helper performs the SciML conversion to a
# `NonlinearProblem` and attaches the sparse Jacobian buffer type
# (`jac_prototype`) that LinearSolve's UMFPACK / KLU factorisations expect.

nlprob = AIBECS.nonlinearproblem(prob)
sol = solve(nlprob, NewtonRaphson(linsolve = UMFPACKFactorization()))

# Or, equivalently, use the AIBECS-recommended default — same Newton step,
# same UMFPACK factorisation:

sol = solve(nlprob, AIBECS.recommended_nlalg())

# ## Sticking with CTKAlg

# `CTKAlg()` continues to work exactly as before — no extra packages
# required:

sol = solve(prob, CTKAlg())

# See the [solvers explanation page](@ref solvers) for guidance on when
# each solver path is appropriate and which algorithm/linear-solver
# combinations are currently verified to work.
