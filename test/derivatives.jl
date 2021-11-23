


# For the conversion to go through ForwardDiff, I need
#Base.convert(FDT::Type{ForwardDiff.Dual{T,V,N}}, p::AIBECS.Parameters) where {T,V,N} = AIBECS.Parameters(convert(Vector{FDT}, vec(p))...)


@testset "Derivatives" begin
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n)
    testp = AIBECS.reconstruct(TestParameters, 2vec(p)) # other value for testing p
    @testset "∇ₓF" begin
        @test ForwardDiff.jacobian(x -> fun(x, p), x) ≈ fun(Val{:jac}, x, p) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> fun(x, p), 2x) ≈ fun(Val{:jac}, 2x, p) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> fun(x, p), -x) ≈ fun(Val{:jac}, -x, p) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> fun(x, testp), x) ≈ fun(Val{:jac}, x, testp) rtol = 1e-14
    end
    @testset "∇ₓf" begin
        @test ForwardDiff.jacobian(x -> [f(x, p)], x) ≈ ∇ₓf(x, p) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> [f(x, testp)], x) ≈ ∇ₓf(x, testp) rtol = 1e-14
        @test ForwardDiff.jacobian(x -> [f(x, p)], 2x) ≈ ∇ₓf(2x, p) rtol = 1e-14
    end
end


