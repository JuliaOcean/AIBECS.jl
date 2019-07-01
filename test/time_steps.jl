
# List of available algorithms for AIBECS time stepping
Implicit_methods = [
                    AIBECS.crank_nicolson_step,
                    AIBECS.crank_nicolson_step!,
                    AIBECS.euler_backward_step,
                    AIBECS.euler_backward_step!,
                   ]

@testset "Implicit time-steps" begin
    nt = length(T_all)
    n = nt * nb
    x₀ = p₀.xgeo * ones(n)
    @testset "$(string(sf))" for sf in Implicit_methods
        δt = ustrip(1.0u"s")
        oldx = copy(x₀)
        x = sf(oldx, p₀, δt, F, ∇ₓF)
        @test x ≈ x₀ + δt * F(x₀, p₀)
    end
end

Explicit_methods = [
                    AIBECS.euler_forward_step,
                    AIBECS.euler_forward_step!,
                   ]

@testset "Explicit time-steps" begin
    nt = length(T_all)
    n = nt * nb
    x₀ = p₀.xgeo * ones(n)
    @testset "$(string(sf))" for sf in Explicit_methods
        δt = ustrip(1.0u"s")
        oldx = copy(x₀)
        x = sf(oldx, p₀, δt, F)
        @test x ≈ x₀ + δt * F(x₀, p₀)
    end
end
