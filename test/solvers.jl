
# List of available algorithms for AIBECS
algs = [CTKAlg]

@testset "Solvers" begin
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n)
    @testset "$(string(alg))" for alg in algs
        @test alg <: DiffEqBase.AbstractSteadyStateAlgorithm
        @testset "$(nx)x, $(np)p" for nx in 1:2, np in 1:2
            testp = AIBECS.reconstruct(typeof(p), np * vec(p))
            prob = SteadyStateProblem(F, ∇ₓF, nx * x, testp)
            s = solve(prob, alg())
            @test s isa SteadyStateSolution
            @test prob isa SteadyStateProblem
            @test norm(s.u) / norm(F(s.u, testp)) > ustrip(upreferred(1e5u"Myr"))
        end
    end
end
