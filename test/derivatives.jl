


# For the conversion to go through ForwardDiff, I need
#Base.convert(FDT::Type{ForwardDiff.Dual{T,V,N}}, p::AIBECS.Parameters) where {T,V,N} = AIBECS.Parameters(convert(Vector{FDT}, vec(p))...)


@testset "Derivatives" begin
    nt = length(T_all)
    n = nt * nb
    x₀ = p₀.xgeo * ones(n)
    @testset "∇ₓF" begin
        @test ForwardDiff.jacobian(x -> F(x, p₀), x₀) ≈ ∇ₓF(x₀, p₀) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> F(x, 2p₀), x₀) ≈ ∇ₓF(x₀, 2p₀) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> F(x, p₀), 2x₀) ≈ ∇ₓF(2x₀, p₀) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> F(x, p₀), -x₀) ≈ ∇ₓF(-x₀, p₀) rtol = 1e-14
    end
    @testset "∇ₚf" begin
        @test ForwardDiff.jacobian(p -> [f(x₀, p)], p₀) ≈ ∇ₚf(x₀, p₀) rtol = 1e-14
        @test ForwardDiff.jacobian(p -> [f(x₀, p)], 2p₀) ≈ ∇ₚf(x₀, 2p₀) rtol = 1e-14
        @test ForwardDiff.jacobian(p -> [f(2x₀, p)], p₀) ≈ ∇ₚf(2x₀, p₀) rtol = 1e-14
    end
    @testset "∇ₓf" begin
        @test ForwardDiff.jacobian(x -> [f(x, p₀)], x₀) ≈ ∇ₓf(x₀, p₀) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> [f(x, 2p₀)], x₀) ≈ ∇ₓf(x₀, 2p₀) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> [f(x, p₀)], 2x₀) ≈ ∇ₓf(2x₀, p₀) rtol = 1e-14
    end
end


