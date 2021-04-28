
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
#=
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
=#


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
sol_SIAMFANL_Newton = let
    F!(dx, x) = F(dx, x, p)
    J!(jac, dx, x) = ∇ₓF(jac, x, p)
    @unpack xgeo = p
    x = xgeo * ones(nb * nt)
    dx = similar(x)
    ftol = ustrip(upreferred(1/1e6Myr))
    nsol(F!, x, dx, fun.jac, J!; rtol=ftol)
end
sol_SIAMFANL_PTC = let
    F!(dx, x) = F(dx, x, p)
    function J2!(jac2, dx, x, dt) # assumes jac2 = jac + diffeqscalar
        for i in 1:length(jac2.ops)-1
            update_coefficients!(jac2.ops[i], x, p, nothing)
        end
        update_coefficients!(jac2.ops[end], nothing, [dt], nothing) # dt part?
        jac2
    end
    function update_dtpart!(A, u, p, t) # here p is just [dt]
        dt = p[1]
        A.diag .= -1 / dt # needs a minus because AIBECS expresses du/dt = F(u) while PTC uses -F(u).
    end
    @unpack xgeo = p
    x = xgeo * ones(nb * nt)
    pdt0=1/norm(F(x,p))
    jac2 = fun.jac + DiffEqArrayOperator(Diagonal(-ones(nb*nt)/pdt0), update_func=update_dtpart!)
    dx = similar(x)
    ftol = ustrip(upreferred(1/1e6Myr))
    ptcsol(F!, x, dx, jac2, J2!; rtol=ftol, pdt0=pdt0, jknowsdt=true, maxit=1000)
end


function fun_SIAMFANL_PTC(dt0)
    F!(dx, x) = F(dx, x, p)
    function J2!(jac2, dx, x, dt) # assumes jac2 = jac + diffeqscalar
        for i in 1:length(jac2.ops)-1
            update_coefficients!(jac2.ops[i], x, p, nothing)
        end
        update_coefficients!(jac2.ops[end], nothing, [dt], nothing) # dt part?
        jac2
    end
    function update_dtpart!(A, u, p, t) # here p is just [dt]
        dt = p[1]
        A.diag .= -1 / dt
    end
    @unpack xgeo = p
    x = xgeo * ones(nb * nt)
    pdt0=dt0/norm(F(x,p))
    jac2 = fun.jac + DiffEqArrayOperator(Diagonal(-ones(nb*nt)/pdt0), update_func=update_dtpart!)
    dx = similar(x)
    ftol = ustrip(upreferred(1/1e6Myr))
    ptcsol(F!, x, dx, jac2, J2!; rtol=ftol, pdt0=pdt0, jknowsdt=true, maxit=10000, keepsolhist=true)
end
