# Hyper parameters
ωs = [1.0, 0.0, 0.0]
ωp = 1e-7
# PO₄ mean and variance of observations fom WOA18
μDIPobs3D, σ²DIPobs3D = WorldOceanAtlasTools.fit_to_grid(grd, "PO₄")
μDIPobs, σ²DIPobs = μDIPobs3D[iwet], σ²DIPobs3D[iwet]
μx = (μDIPobs, missing, missing)
σ²x = (σ²DIPobs, missing, missing)
δconvert(x) = @. exp(x)
cs = (δconvert, identity, identity)
# generate mismatch functions
f, ∇ₓf = f_and_∇ₓf(ωs, μx, σ²x, v, ωp, TestParameters)

#===========================================
Tests
===========================================#

# TODO add tests for other types of mismatch functions!

@testset "Mismatch / cost functions" begin
    nt = length(T_all)
    nx = nt * nb
    np = length(p)
    x = p.xgeo * ones(nx)
    λ = p2λ(p)
    fλ = f(x, λ)
    ∇ₓfλ = ∇ₓf(x, λ)
    @testset "Mismatch function" begin
        @test fλ isa Float64
        @test fλ > 0
        @test fλ ≈ f(x, p)
    end
    @testset "Derivative w.r.t the state, x" begin
        @test size(∇ₓfλ) == (1, nx)
        @test ∇ₓfλ ≈ ∇ₓf(x, p)
    end
end
