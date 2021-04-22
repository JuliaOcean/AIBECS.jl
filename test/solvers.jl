
# List of available algorithms for AIBECS
#algs = [CTKAlg]


#@testset "Solvers" begin
#    n = nt * nb
#    @unpack xgeo = p
#    x = xgeo * ones(n)
#    @testset "$(string(alg))" for alg in algs
#        @test alg <: DiffEqBase.AbstractSteadyStateAlgorithm
#        @testset "$(nx)x, $(np)p" for nx in 1:2, np in 1:2
#            testp = AIBECS.reconstruct(TestParameters, np * vec(p))
#            oopprob = SteadyStateProblem(fun, nx * x, testp)
#            s = solve(oopprob, alg())
#            @test s isa SteadyStateSolution
#            @test oopprob isa SteadyStateProblem
#            @test norm(s.u) / norm(fun.f(s.u, testp, 0)) > ustrip(upreferred(1e5u"Myr"))
#
#            iipprob = SteadyStateProblem(fun!, nx * x, testp)
#            s = solve(iipprob, alg())
#            @test s isa SteadyStateSolution
#            @test iipprob isa SteadyStateProblem
#            @test norm(s.u) / norm(fun!.f(s.u, testp, 0)) > ustrip(upreferred(1e5u"Myr"))
#        end
#    end
#end




#===========================================
my CTKAlg
===========================================#
sol_CTKAlg = let
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
            @test norm(s.u) / norm(fun!.f(s.u, testp, 0)) > ustrip(upreferred(1e5u"Myr"))
        end
    end
end


#===========================================
CTKAlg
===========================================#
sol_CTKAlg = let
    f(x) = F(x, p)
    j(x) = ∇ₓF(x, p)
    @unpack xgeo = p
    x = xgeo * ones(nb * nt)
    τstop = ustrip(upreferred(1e6Myr))
    AIBECS.NewtonChordShamanskii(f, j, norm, x, τstop)
end


#===========================================
NLsolve
===========================================#
sol_NLsolve = let
    f(x) = F(x, p)
    j(x) = ∇ₓF(x, p)
    @unpack xgeo = p
    x = xgeo * ones(nb * nt)
    method = :newton
    ftol = ustrip(upreferred(1/1e6Myr))
    nlsolve(f, j, x; method, ftol)
end


#===========================================
SIAMFANLEquations
===========================================#
sol_NLSolvers = let
    F!(dx, x) = F(dx, x, p)
    J!(jac, dx, x) = ∇ₓF(jac, x, p)
    @unpack xgeo = p
    x = xgeo * ones(nb * nt)
    dx = similar(x)
    ftol = ustrip(upreferred(1/1e6Myr))
    nsol(F!, x, dx, fun.jac, J!; rtol=ftol)
end
sol_NLSolvers_PTC = let
    F!(dx, x) = F(dx, x, p)
    J!(jac, dx, x) = ∇ₓF(jac, x, p)
    @unpack xgeo = p
    x = xgeo * ones(nb * nt)
    dx = similar(x)
    ftol = ustrip(upreferred(1/1e6Myr))
    ptcsol(F!, x, dx, fun.jac, J!; rtol=ftol)
end

+(::DiffEqOpCombination{
                              Float64, 
                              Tuple{
                                    DiffEqArrayOp{Float64, SparseMatrixCSC{Float64, Int64}, AIBECS.var"#update_∇G#87"{AIBECS.var"#G#80"{Tuple{typeof(sms_DIP), typeof(sms_DOP), typeof(sms_POP)}, Int64, AIBECS.var"#tracers#77"{Int64, Int64}, AIBECS.var"#tracer#78"{Int64, Int64}}, Vector{Int64}, SparseMatrixCSC{Float64, Int64}, Nothing}},
                                    DiffEqScaledOp{Float64, typeof(SciMLBase.DEFAULT_UPDATE_FUNC), AIBECS.BlockDiagDiffEqOp{Float64, Tuple{DiffEqArrayOp{Float64, SparseMatrixCSC{Float64, Int64}, typeof(SciMLBase.DEFAULT_UPDATE_FUNC)}, DiffEqArrayOp{Float64, SparseMatrixCSC{Float64, Int64}, typeof(SciMLBase.DEFAULT_UPDATE_FUNC)}, DiffEqArrayOp{Float64, SparseMatrixCSC{Float64, Int64}, AIBECS.var"#update_func#11"{typeof(w), Vector{Float64}, AIBECS.var"#10#17", Vector{Float64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Int64}, Vector{Float64}, Vector{Int64}}}}}}
                                   },
                              Vector{Float64}
                             },
  ::UniformScaling{Float64})


