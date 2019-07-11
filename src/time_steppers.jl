# Some time-steppers for AIBECS simulations

# TODO:
# - increment other methods
#     - IMEX
#     - other bashford something methods
# - Try to integrate into DiffEq implementation instead

"""
    euler_forward_step(x, p, δt, F)

Returns the Euler-forward increment `x + δt * F(x,p)` (out of place).
"""
euler_forward_step(x, p, δt, F) = x + δt * F(x, p)
"""
    euler_forward_step!(x, p, δt, F)

In place version of [`euler_forward_step`](@ref).
"""
euler_forward_step!(x, p, δt, F) = x .= euler_forward_step(x, p, δt, F)

"""
    euler_backward_step(x, p, δt, F, ∇ₓF)

Returns the Euler-backward increment (out of place).
Solves `H(y) = 0` where `H(y) = F(y,p) + (x - y) / δt`.
"""
function euler_backward_step(x, p, δt, F, ∇ₓF)
    H(y,p) = F(y,p) - y / δt + x / δt
    ∇ₓH(y,p) = ∇ₓF(y,p) - I / δt
    prob = SteadyStateProblem(H, ∇ₓH, x, p)
    return solve(prob, CTKAlg()).u
end
"""
    euler_backward_step!(x, p, δt, F, ∇ₓF)

In place version of [`euler_backward_step`](@ref).
"""
function euler_backward_step!(x, p, δt, F, ∇ₓF)
    H(y,p) = F(y,p) - y / δt + x / δt
    ∇ₓH(y,p) = ∇ₓF(y,p) - I / δt
    prob = SteadyStateProblem(H, ∇ₓH, x, p)
    x .= solve(prob, CTKAlg()).u
    return x
end

"""
    crank_nicolson_step(x, p, δt, F, ∇ₓF)

Returns the crank nicolson increment (out of place).
Solves `H(y) = 0` where `H(y) = (F(y,p) + F(x,p)) / 2 + (x - y) / δt`.
"""
function crank_nicolson_step(x, p, δt, F, ∇ₓF)
    H(y,p) = F(y,p) / 2 - y / δt + F(x,p) / 2 + x / δt
    ∇ₓH(y,p) = ∇ₓF(y,p) / 2 - I / δt
    prob = SteadyStateProblem(H, ∇ₓH, x, p)
    return solve(prob, CTKAlg()).u
end
"""
    crank_nicolson_step!(x, p, δt, F, ∇ₓF)

In place version of [`crank_nicolson_step`](@ref).

For a tracer equation of the form

`∂x/∂t = F(x,p)`

The Crank-Nicolson scheme is given by

`(xᵢ₊₁ - xᵢ) / δt = (F(xᵢ₊₁,p) + F(xᵢ,p)) / 2`

which requires solving `H(y) = 0` where `H` is defined by

`H(y) = (F(y,p) + F(xᵢ,p)) / 2 + (xᵢ - y) / δt`

where `xᵢ` is prescribed.
`H(y) = 0` is thus defined as a steady-state problem and solved via Newton's method.
"""
function crank_nicolson_step!(x, p, δt, F, ∇ₓF)
    H(y,p) = F(y,p) / 2 - y / δt + F(x,p) / 2 + x / δt
    ∇ₓH(y,p) = ∇ₓF(y,p) / 2 - I / δt
    prob = SteadyStateProblem(H, ∇ₓH, x, p)
    x .= solve(prob, CTKAlg()).u
    return x
end

"""
    crank_nicolson_leapfrog_step(xᵢ, xᵢ₋₁, p, δt, T, ∇ₓL, NL)

Returns the Crank-Nicolson-leapfrog-step increment (out of place).
See [`crank_nicolson_leapfrog_step_A⁺_and_A⁻`](@ref) for more details and a more efficient
computation.
"""
function crank_nicolson_leapfrog_step(xᵢ, xᵢ₋₁, p, δt, T, ∇ₓL, NL)
    return (I / 2δt + (T(p) - ∇ₓL(p)) / 2) \ (NL(xᵢ,p) + (I / 2δt - (T(p) - ∇ₓL(p)) / 2) * xᵢ₋₁)
end
"""
    crank_nicolson_leapfrog_step(xᵢ, xᵢ₋₁, p, A⁺::Factorization, A⁻, G)

Returns the Crank-Nicolson-leapfrog-step increment (out of place) using stored factors.
See [`crank_nicolson_leapfrog_step_A⁺_and_A⁻`](@ref) for more information.
"""
function crank_nicolson_leapfrog_step(xᵢ, xᵢ₋₁, p, A⁺::Factorization, A⁻, NL)
    return A⁺ \ (NL(xᵢ,p) + A⁻ * xᵢ₋₁)
end
"""
    crank_nicolson_leapfrog_step_A⁺_and_A⁻(p, δt, T, ∇ₓL)

Returns `A⁺` and `A⁻` for the Crank-Nicolson-leapfrog-step increment (out of place).
`A⁺` is a `Factorization` object and `A⁻` is a matrix.

For a tracer equation of the form

`∂x/∂t = F(x,p) = NL(x,p) + [∇ₓL(p) - T(p)] x`

where `F(x,p)` has been split into 3 terms:
- the linear transport, `T(p) x`,
- a non linear source/sink term, `NL(x,p)`,
- and a linear source/sink term, `L(x,p) = ∇ₓL(p) x` (note the signs of `∇ₓL` and `T`).
The Crank-Nicolson-leapfrog scheme uses this nonlinear/linear split and is given by

`(xᵢ₊₁ - xᵢ₋₁) / 2δt = NL(xᵢ,p) + [∇ₓL(p) - T(p)] (xᵢ₊₁ + xᵢ₋₁) / 2`

or, with the `xᵢ₊₁` terms on the left,

`A⁺ xᵢ₊₁ = NL(xᵢ,p) + A⁻ xᵢ₋₁`

with `2A⁺ = I / δt + T(p) - ∇ₓL(p)` and `2A⁻ = I / δt - T(p) + ∇ₓL(p)`,
which only depend on `p` (not `x`).
Thus, the Crank-Nicolson-leapfrog steps only needs to factorize `A⁺` once, making it much faster than the standard Crank-Nicolson scheme, justifying this function.
"""
function crank_nicolson_leapfrog_step_A⁺_and_A⁻(p, δt, T, ∇ₓL)
    return factorize(I / 2δt + (T(p) - ∇ₓL(p)) / 2), I / 2δt - (T(p) - ∇ₓL(p)) / 2
end
