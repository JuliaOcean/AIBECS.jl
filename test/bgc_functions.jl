#===========================================
Transport matrices
===========================================#
T_DIP(p) = T_Circulation
T_DOP(p) = T_Circulation
S₀ = buildPFD(ones(nb), DIV, Iabove)
S′ = buildPFD(ztop, DIV, Iabove)
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
    τgeo, xgeo = p.τgeo, p.xgeo
    return @. (xgeo - x) / τgeo
end
# Uptake of phosphate (DIP)
relu(x) = (x .≥ 0) .* x
function uptake(DIP, p)
    τ, k, z₀ = p.τ, p.k, p.z₀
    DIP⁺ = relu(DIP)
    return @. 1 / τ * DIP⁺^2 / (DIP⁺ + k) * (z ≤ z₀)
end
# Remineralization DOP into DIP
function remineralization(DOP, p)
    κDOP = p.κDOP
    return @. κDOP * DOP
end
# Dissolution of POP into DOP
function dissolution(POP, p)
    κPOP = p.κPOP
    return @. κPOP * POP
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
In-place sources minus sinks
===========================================#

function G_DIP!(dx, DIP, DOP, POP, p)
    τgeo, xgeo, τ, k, z₀, κDOP = p.τgeo, p.xgeo, p.τ, p.k, p.z₀, p.κDOP
    dx .= @. (xgeo - DIP) / τgeo - (DIP ≥ 0) / τ * DIP^2 / (DIP + k) * (z ≤ z₀) + κDOP * DOP
end
function G_DOP!(dx, DIP, DOP, POP, p)
    τgeo, xgeo, τ, k, z₀, κDOP, κPOP, σ = p.τgeo, p.xgeo, p.τ, p.k, p.z₀, p.κDOP, p.κPOP, p.σ
    dx = @. σ * (DIP ≥ 0) / τ * DIP^2 / (DIP + k) * (z ≤ z₀) - κDOP * DOP + κPOP * POP
end
function G_POP!(dx, DIP, DOP, POP, p)
    τgeo, xgeo, τ, k, z₀, κDOP, κPOP, σ = p.τgeo, p.xgeo, p.τ, p.k, p.z₀, p.κDOP, p.κPOP, p.σ
    dx = @. (1 - σ) * (DIP ≥ 0) / τ * DIP^2 / (DIP + k) * (z ≤ z₀) - κPOP * POP
end
Gs = (G_DIP!, G_DOP!, G_POP!)

#===========================================
AIBECS F and ∇ₓF
===========================================#
F, ∇ₓF = state_function_and_Jacobian(T_all, sms_all, nb)
F!, ∇ₓF! = inplace_state_function_and_Jacobian(T_all, Gs, nb)

#===========================================
AIBECS split operators for CNLF
===========================================#

NL_DIP(DIP, DOP, POP, p) = -uptake(DIP, p) + geores(0*DIP, p)
L_DIP(DIP, DOP, POP, p) = remineralization(DOP, p) + geores(DIP, p) - geores(0*DIP, p)

NL_DOP(DIP, DOP, POP, p) = p.σ * uptake(DIP, p)
L_DOP(DIP, DOP, POP, p) = -remineralization(DOP, p) + dissolution(POP, p)

NL_POP(DIP, DOP, POP, p) = (1 - p.σ) * uptake(DIP, p)
L_POP(DIP, DOP, POP, p) = -dissolution(POP, p)

NL_all = (NL_DIP, NL_DOP, NL_POP)
L_all = (L_DIP, L_DOP, L_POP)

F2, L, NL, ∇ₓF2, ∇ₓL, ∇ₓNL, T = split_state_function_and_Jacobian(T_all, L_all, NL_all, nb)

#===========================================
Tests
===========================================#

@testset "Biogeochemical functions" begin
    @testset "Particle Flux Divergence is conservative" begin
        @test norm(v) / norm(T_POP(p₀)' * v) > ustrip(upreferred(1u"Myr"))
    end
    nt = length(T_all)
    n = nt * nb
    x = p₀.xgeo * ones(n) .* exp.(cos.(collect(1:n))/10)
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


@testset "In-place functions" begin
    nt = length(T_all)
    n = nt * nb
    x = p₀.xgeo * ones(n) .* exp.(cos.(collect(1:n))/10)
    dx = similar(x)
    @testset "F! ≈ F" begin
        F₀ = F!(dx, x, p₀)
        @test F₀ ≈ F(x,p₀) rtol=1e-10
    end
    @testset "∇ₓF! ≈ ∇ₓF" begin
        @test ∇ₓF(x,p₀) ≈ ∇ₓF!(x,p₀) rtol=1e-10
    end
end

@testset "Linear/nonlinear split functions" begin
    nt = length(T_all)
    n = nt * nb
    x = p₀.xgeo * ones(n) .* exp.(cos.(collect(1:n))/10)
    @testset "split F ≈ F" begin
        @test F2(x,p₀) ≈ F(x,p₀) rtol=1e-10
    end
    @testset "split F ≈ L + NL" begin
        @test F2(x,p₀) ≈ L(x,p₀) + NL(x,p₀) - T(p₀) * x rtol=1e-10
    end
    @testset "split ∇ₓF ≈ ∇ₓF" begin
        @test ∇ₓF(x,p₀) ≈ ∇ₓF2(x,p₀) rtol=1e-10
    end
    @testset "split ∇ₓF ≈ T + ∇ₓNL" begin
        @test ∇ₓF2(x,p₀) ≈ ∇ₓL(p₀) + ∇ₓNL(x,p₀) - T(p₀) rtol=1e-10
    end
end
