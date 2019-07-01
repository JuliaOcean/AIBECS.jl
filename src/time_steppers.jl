# Some time-steppers for AIBECS simulations

# TODO:
# - increment other methods
#     - IMEX_step(x, δt, F, ) = (I + δt * M) \ (R + δt * s)
#     - crank_nicolson_leapfrog_step(xᵢ, xᵢ₋₁, δt, L, G) = xᵢ₋₁ + L \ G(xᵢ)
# - Try to integrate into DiffEq implementation instead

"""
    euler_forward_step(x, p, δt, F)

Returns the Euler-forward increment `x + δt * F(x,p)` (out of place).
"""
euler_forward_step(x, p, δt, F) = x + δt * F(x, p)
"""
    euler_forward_step!(x, p, δt, F)

Returns the Euler-forward increment `x + δt * F(x,p)` (updates `x` in place).
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

Returns the Euler-backward increment (updates `x` in place).
Solves `H(y) = 0` where `H(y) = F(y,p) + (x - y) / δt`.
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

Returns the crank nicolson increment (updates `x` in place).
Solves `H(y) = 0` where `H(y) = (F(y,p) + F(x,p)) / 2 + (x - y) / δt`.
"""
function crank_nicolson_step!(x, p, δt, F, ∇ₓF)
    H(y,p) = F(y,p) / 2 - y / δt + F(x,p) / 2 + x / δt
    ∇ₓH(y,p) = ∇ₓF(y,p) / 2 - I / δt
    prob = SteadyStateProblem(H, ∇ₓH, x, p)
    x .= solve(prob, CTKAlg()).u
    return x
end

