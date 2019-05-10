# Hyper parameters
ωs = [1.0, 0.0]
ωp = 1e-4
# PO₄ mean and variance of observations fom WOA18
WOA = WorldOceanAtlasTools
μDIPobs3D, σ²DIPobs3D = WOA.fit_to_grid(grd, 2018, "phosphate", "annual", "1°", "an")
μDIPobs, σ²DIPobs = μDIPobs3D[iwet], σ²DIPobs3D[iwet]
const μx = (μDIPobs, missing)
const σ²x = (σ²DIPobs, missing)
# generate mismatch functions
f = generate_objective(ωs, μx, σ²x, v, ωp, mean_pobs, variance_pobs)
∇ₓf = generate_∇ₓobjective(ωs, μx, σ²x, v, ωp, mean_pobs, variance_pobs)
∇ₚf = generate_∇ₚobjective(ωs, μx, σ²x, v, ωp, mean_pobs, variance_pobs)

#===========================================
Tests
===========================================#

@testset "Mismatch / cost functions" begin
    x₀ = p₀.DIPgeo * kron([1, 0.1], ones(nb))
    nt = length(T_all)
    n = nt * nb
    f₀ = f(x₀, p₀)
    ∇ₓf₀ = ∇ₓf(x₀, p₀)
    ∇ₚf₀ = ∇ₚf(x₀, p₀)
    @testset "Mismatch function" begin
        @test f₀ isa Float64
        @test f₀ > 0
    end
    @testset "Derivative w.r.t the state, x" begin
        @test size(∇ₓf₀) == (1, n)
    end
    @testset "Derivative w.r.t the parameters, p" begin
        @test size(∇ₚf₀) == (1, m)
    end
end
