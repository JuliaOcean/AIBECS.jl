


# For the conversion to go through ForwardDiff, I need
#Base.convert(FDT::Type{ForwardDiff.Dual{T,V,N}}, p::AIBECS.Parameters) where {T,V,N} = AIBECS.Parameters(convert(Vector{FDT}, vec(p))...)


@testset "Derivatives" begin
    x₀ = p₀.DIPgeo * kron([1, 0.1], ones(nb))
    @testset "∇ₓF" begin
        @test ForwardDiff.jacobian(x -> F(x, p₀), x₀) ≈ ∇ₓF(x₀, p₀) rtol = 1e-12
    end
    @testset "∇ₚf" begin
        @test ForwardDiff.jacobian(p -> [f(x₀, p)], p₀) ≈ ∇ₚf(x₀, p₀) rtol = 1e-12
    end
    @testset "∇ₓf" begin
        @test ForwardDiff.jacobian(x -> [f(x, p₀)], x₀) ≈ ∇ₓf(x₀, p₀) rtol = 1e-12
    end
end


