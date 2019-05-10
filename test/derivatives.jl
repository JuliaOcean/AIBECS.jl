# For the conversion to go through ForwardDiff, I need
Base.convert(FDT::Type{ForwardDiff.Dual{T,V,N}}, p::AIBECS.Parameters) where {T,V,N} = AIBECS.Parameters(convert(Vector{FDT}, vec(p))...)

@test ForwardDiff.jacobian(v -> F(x₀, AIBECS.opt_para(v)), optvec(p₀)) ≈ ∇ₚF(x₀, p₀)
@test ForwardDiff.jacobian(x -> F(x, p₀), x₀) ≈ ∇ₓF(x₀, p₀)
@test ForwardDiff.jacobian(p -> f(x₀, p), p₀) ≈ ∇ₚf(x₀, p₀)
@test ForwardDiff.jacobian(x -> f(x, p₀), x₀) ≈ ∇ₓf(x₀, p₀)


