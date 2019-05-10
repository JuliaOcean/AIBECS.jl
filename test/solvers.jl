
# List of available algorithms for AIBECS
algs = [CTKAlg]

@testset "Solvers" begin
    @testset "$string(alg)" for alg in algs
        x₀ = p₀.DIPgeo * kron([1, 0.1], ones(nb))
        prob = SteadyStateProblem(F, ∇ₓF, x₀, p₀)
        s = solve(prob, alg())
        @test alg <: DiffEqBase.AbstractSteadyStateAlgorithm
        @test s isa SteadyStateSolution
        @test prob isa SteadyStateProblem
        @test norm(s.u) / norm(F(s.u, p₀)) > ustrip(upreferred(1u"Myr"))
    end
end
