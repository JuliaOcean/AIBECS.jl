
# List of available algorithms for AIBECS
algs = [CTKAlg]

@testset "Solvers" begin
    @testset "$string(alg)" for alg in algs
        x₀ = p₀.DIPgeo * kron([1, 0.1], ones(nb))
        @test alg <: DiffEqBase.AbstractSteadyStateAlgorithm
        @testset "$(nx₀)x₀, $(np₀)p₀" for nx₀ in 1:2, np₀ in 1:2
            prob = SteadyStateProblem(F, ∇ₓF, nx₀ * x₀, np₀ * p₀)
            s = solve(prob, alg())
            @test s isa SteadyStateSolution
            @test prob isa SteadyStateProblem
            @test norm(s.u) / norm(F(s.u, np₀ * p₀)) > ustrip(upreferred(1e5u"Myr"))
        end
    end
end
