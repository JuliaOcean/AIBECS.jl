
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
            testp = AIBECS.reconstruct(TestParameters, np * vec(p))
            oopprob = SteadyStateProblem(fun, nx * x, testp)
            s = solve(oopprob, alg())
            @test s isa SteadyStateSolution
            @test oopprob isa SteadyStateProblem
            @test norm(s.u) / norm(fun.f(s.u, testp, 0)) > ustrip(upreferred(1e5u"Myr"))

            iipprob = SteadyStateProblem(fun!, nx * x, testp)
            s = solve(iipprob, alg())
            @test s isa SteadyStateSolution
            @test iipprob isa SteadyStateProblem
            ds = copy(s.u)
            @test norm(s.u) / norm(fun!.f(ds, s.u, testp, 0)) > ustrip(upreferred(1e5u"Myr"))
        end
    end
end
