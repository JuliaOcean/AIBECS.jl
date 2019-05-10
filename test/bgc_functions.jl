#===========================================
Transport matrices
===========================================#
T_DIP(p) = T_Circulation
const S₀ = buildPFD(ones(nb), DIV, Iabove)
const S′ = buildPFD(ztop, DIV, Iabove)
function T_POP(p)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S₀ + w′ * S′
end
T_all = (T_DIP, T_POP)

#===========================================
Sources minus sinks
===========================================#
# Geological Restoring
function geores(DIP, p)
    τg, DIPgeo = p.τg, p.DIPgeo
    return (DIPgeo .- DIP) ./ τg
end
# Uptake of phosphate (DIP)
relu(x) = (x .≥ 0) .* x
function uptake(DIP, p)
    Umax, ku, z₀ = p.Umax, p.ku, p.z₀
    DIP⁺ = relu(DIP)
    return Umax * DIP⁺ ./ (DIP⁺ .+ ku) .* (z .≤ z₀)
end
# Remineralization of particulate organic phosphorus (POP)
function remineralization(POP, p)
    κ = p.κ
    return κ * POP
end
# Add them up into sms functions (Sources Minus Sinks)
sms_DIP(DIP, POP, p) = geores(DIP, p) .- uptake(DIP, p) .+ remineralization(POP, p)
sms_POP(DIP, POP, p) = uptake(DIP, p) .- remineralization(POP, p)
sms_all = (sms_DIP, sms_POP) # bundles all the source-sink functions in a tuple

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
    x = p₀.DIPgeo * kron([1, 0.1], ones(nb))
    nt = length(T_all)
    n = nt * nb
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
