
# List of available algorithms for AIBECS time stepping

@testset "Time steppers" begin
    δt = ustrip(1.0u"yr" |> u"s")
    Nt = 10000 # number of time steps taken
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x₀ = xgeo * ones(n)

    Implicit_methods = [
                        AIBECS.crank_nicolson_step,
                        AIBECS.crank_nicolson_step!,
                        AIBECS.euler_backward_step,
                        AIBECS.euler_backward_step!,
                       ]
    @testset "Implicit time steps" begin
        @testset "$(string(sf))" for sf in Implicit_methods
            x = copy(x₀)
            for i in 1:Nt
                x .= sf(x, p, δt, F, ∇ₓF)
            end
            @test norm(x) / norm(F(x,p)) > ustrip(upreferred(1u"kyr"))
        end
    end

    Explicit_methods = [
                        AIBECS.euler_forward_step,
                        AIBECS.euler_forward_step!,
                       ]
    @testset "Explicit time steps" begin
        @testset "$(string(sf))" for sf in Explicit_methods
            x = copy(x₀)
            for i in 1:100Nt              # more smaller steps
                x .= sf(x, p, 0.01δt, F) # for explicit scheme
            end
            @test norm(x) / norm(F(x,p)) > ustrip(upreferred(1u"kyr"))
        end
    end

    @testset "Crank-Nicolson leapfrog" begin
        x, xᵢ, xᵢ₋₁ = copy(x₀), copy(x₀), copy(x₀)
        for i in 1:Nt
            x .= AIBECS.crank_nicolson_leapfrog_step(xᵢ, xᵢ₋₁, p, δt, T, ∇ₓL, NL)
            xᵢ₋₁ .= xᵢ
            xᵢ .= x
        end
        @test norm(x) / norm(F(x,p)) > ustrip(upreferred(1u"kyr"))
        A⁺, A⁻ = AIBECS.crank_nicolson_leapfrog_step_A⁺_and_A⁻(p, δt, T, ∇ₓL)
        x, xᵢ, xᵢ₋₁ = copy(x₀), copy(x₀), copy(x₀)
        for i in 1:Nt
            x .= AIBECS.crank_nicolson_leapfrog_step(xᵢ, xᵢ₋₁, p, A⁺, A⁻, NL)
            xᵢ₋₁ .= xᵢ
            xᵢ .= x
        end
        @test norm(x) / norm(F(x,p)) > ustrip(upreferred(1u"kyr"))
    end

end
