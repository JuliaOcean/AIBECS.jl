# Hyper parameters
ωs = [1.0, 0.0, 0.0]
ωp = 1e-7
# PO₄ mean and variance of observations fom WOA18
WOA = WorldOceanAtlasTools
μDIPobs3D, σ²DIPobs3D = WOA.fit_to_grid(grid, 2018, "phosphate", "annual", "1°", "an")
μDIPobs, σ²DIPobs = μDIPobs3D[iwet], σ²DIPobs3D[iwet]
μx = (μDIPobs, missing, missing)
σ²x = (σ²DIPobs, missing, missing)
# generate mismatch functions
f, ∇ₓf, ∇ₚf = generate_objective_and_derivatives(ωs, μx, σ²x, v, ωp)

#===========================================
Tests
===========================================#

@testset "Mismatch / cost functions" begin
    nt = length(T_all)
    n = nt * nb
    m = length(p)
    x = p.xgeo * ones(n)
    fval = f(x, p)
    ∇ₓfval = ∇ₓf(x, p)
    ∇ₚfval = ∇ₚf(x, p)
    @testset "Mismatch function" begin
        @test fval isa Float64
        @test fval > 0
    end
    @testset "Derivative w.r.t the state, x" begin
        @test size(∇ₓfval) == (1, n)
    end
    @testset "Derivative w.r.t the parameters, p" begin
        @test size(∇ₚfval) == (1, m)
    end
end
