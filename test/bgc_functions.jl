#===========================================
Transport matrices
===========================================#
T_DIP(p) = T_Circulation
T_DOP(p) = T_Circulation
const S₀ = buildPFD(ones(nb), DIV, Iabove)
const S′ = buildPFD(ztop, DIV, Iabove)
function T_POP(p)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S₀ + w′ * S′
end
T_all = (T_DIP, T_DOP, T_POP)

#===========================================
Sources minus sinks
===========================================#
# Geological Restoring
function geores(x, p)
    τg, xgeo = p.τg, p.xgeo
    return (xgeo .- x) / τg
end
# Uptake of phosphate (DIP)
soft_relu(x, p) = 0.5 * (tanh.(x / p.α) .+ 1) .* x
function uptake(DIP, p)
    Umax, ku, z₀ = p.Umax, p.ku, p.z₀
    DIP⁺ = soft_relu(DIP, p)
    return Umax * DIP⁺ ./ (DIP⁺ .+ ku) .* (z .≤ z₀)
end
# Remineralization DOP into DIP
function remineralization(DOP, p)
    κDOP = p.κDOP
    return κDOP * DOP
end
# Dissolution of POP into DOP
function dissolution(POP, p)
    κPOP = p.κPOP
    return κPOP * POP
end
# Add them up into sms functions (Sources Minus Sinks)
function sms_DIP(DIP, DOP, POP, p)
    return -uptake(DIP, p) + remineralization(DOP, p) + geores(DIP, p)
end
function sms_DOP(DIP, DOP, POP, p)
    σ = p.σ
    return σ * uptake(DIP, p) - remineralization(DOP, p) + dissolution(POP, p)
end
function sms_POP(DIP, DOP, POP, p)
    σ = p.σ
    return (1 - σ) * uptake(DIP, p) - dissolution(POP, p)
end
sms_all = (sms_DIP, sms_DOP, sms_POP) # bundles all the source-sink functions in a tuple

#===========================================
AIBECS F and ∇ₓF
===========================================#
F, ∇ₓF = state_function_and_Jacobian(T_all, sms_all, nb)

#===========================================
Tests
===========================================#

@testset "Biogeochemical functions" begin
    @testset "Particle Flux Divergence is conservative" begin
        @test norm(v) / norm(T_POP(p₀)' * v) > ustrip(upreferred(1u"Myr"))
    end
    nt = length(T_all)
    n = nt * nb
    x = p₀.xgeo * ones(n)
    @testset "State function" begin
        F₀ = F(x, p₀)
        v_all = repeat(v, nt)
        @test F₀ isa Vector{Float64}
        @test size(F₀) == (n,)
        @test norm(v_all .* x) / norm(v_all' * F₀) > ustrip(upreferred(1u"Myr"))
    end
    @testset "Jacobian of the state function" begin
        ∇ₓF₀ = ∇ₓF(x, p₀)
        @test ∇ₓF₀ isa SparseMatrixCSC
        @test size(∇ₓF₀) == (n,n)
    end
end
