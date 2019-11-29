


# For the conversion to go through ForwardDiff, I need
#Base.convert(FDT::Type{ForwardDiff.Dual{T,V,N}}, p::AIBECS.Parameters) where {T,V,N} = AIBECS.Parameters(convert(Vector{FDT}, vec(p))...)


@testset "Derivatives" begin
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n)
    testp = reconstruct(typeof(p), 2vec(p)) # other value for testing p
    @testset "∇ₓF" begin
        @test ForwardDiff.jacobian(x -> F(x, p), x) ≈ ∇ₓF(x, p) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> F(x, p), 2x) ≈ ∇ₓF(2x, p) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> F(x, p), -x) ≈ ∇ₓF(-x, p) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> F(x, testp), x) ≈ ∇ₓF(x, testp) rtol = 1e-14
    end
    @testset "∇ₚf" begin
        @test ForwardDiff.jacobian(p -> [f(x, p)], p) ≈ ∇ₚf(x, p) rtol = 1e-14
        @test ForwardDiff.jacobian(p -> [f(x, p)], testp) ≈ ∇ₚf(x, testp) rtol = 1e-14
        @test ForwardDiff.jacobian(p -> [f(2x, p)], p) ≈ ∇ₚf(2x, p) rtol = 1e-14
    end
    @testset "∇ₓf" begin
        @test ForwardDiff.jacobian(x -> [f(x, p)], x) ≈ ∇ₓf(x, p) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> [f(x, testp)], x) ≈ ∇ₓf(x, testp) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> [f(x, p)], 2x) ≈ ∇ₓf(2x, p) rtol = 1e-14
    end
end


