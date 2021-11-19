#===========================================
Transport matrices
===========================================#
T_D(p) = T_Circulation
S₀ = PFDO(ones(nb), DIV, Iabove)
S′ = PFDO(ztop, DIV, Iabove)
function T_POP(p)
    @unpack w₀, w′ = p
    return w₀ * S₀ + w′ * S′
end
T_all = (T_D, T_D, T_POP)

#===========================================
Sources minus sinks
===========================================#
# Geological Restoring
function geores(x, p)
    @unpack τgeo, xgeo = p
    return @. (xgeo - x) / τgeo
end
# Uptake of phosphate (DIP)
relu(x) = (x .≥ 0) .* x
function uptake(DIP, p)
    @unpack τDIP, k = p
    DIP⁺ = relu(DIP)
    return @. 1 / τDIP * DIP⁺^2 / (DIP⁺ + k) * (ztop == 0)
end
# Remineralization DOP into DIP
function remineralization(DOP, p)
    @unpack τDOP = p
    return @. DOP / τDOP
end
# Dissolution of POP into DOP
function dissolution(POP, p)
    @unpack τPOP = p
    return @. POP / τPOP
end
# Add them up into sms functions (Sources Minus Sinks)
function sms_DIP(DIP, DOP, POP, p)
    return -uptake(DIP, p) + remineralization(DOP, p) + geores(DIP, p)
end
function sms_DOP(DIP, DOP, POP, p)
    @unpack σ = p
    return σ * uptake(DIP, p) - remineralization(DOP, p) + dissolution(POP, p)
end
function sms_POP(DIP, DOP, POP, p)
    @unpack σ = p
    return (1 - σ) * uptake(DIP, p) - dissolution(POP, p)
end
sms_all = (sms_DIP, sms_DOP, sms_POP) # bundles all the source-sink functions in a tuple

#===========================================
In-place sources minus sinks
===========================================#

function G_DIP!(dx, DIP, DOP, POP, p)
    @unpack τgeo, xgeo, τDIP, k, τDOP = p
    @. dx = (xgeo - DIP) / τgeo - (DIP ≥ 0) / τDIP * DIP^2 / (DIP + k) * (ztop == 0) + DOP / τDOP
end
function G_DOP!(dx, DIP, DOP, POP, p)
    @unpack τgeo, xgeo, τDIP, k, τDOP, τPOP, σ = p
    @. dx = σ * (DIP ≥ 0) / τDIP * DIP^2 / (DIP + k) * (ztop == 0) - DOP / τDOP + POP / τPOP
end
function G_POP!(dx, DIP, DOP, POP, p)
    @unpack τgeo, xgeo, τDIP, k, τDOP, τPOP, σ = p
    @. dx = (1 - σ) * (DIP ≥ 0) / τDIP * DIP^2 / (DIP + k) * (ztop == 0) - POP / τPOP
end
Gs = (G_DIP!, G_DOP!, G_POP!)


#===========================================
AIBECS F and ∇ₓF
===========================================#
fun = AIBECSFunction(T_all, sms_all, nb)
fun! = AIBECSFunction((T_D, T_POP), Gs, nb, 3, [1, 1, 2])
F, ∇ₓF = F_and_∇ₓF(T_all, sms_all, nb) # deprecated
F!, ∇ₓF! = F_and_∇ₓF((T_D, T_POP), Gs, nb, 3, [1, 1, 2]) # deprecated

#===========================================
AIBECS split operators for CNLF
===========================================#

NL_DIP(DIP, DOP, POP, p) = -uptake(DIP, p) + geores(0 * DIP, p)
L_DIP(DIP, DOP, POP, p) = remineralization(DOP, p) + geores(DIP, p) - geores(0 * DIP, p)

function NL_DOP(DIP, DOP, POP, p)
    @unpack σ = p
    return σ * uptake(DIP, p)
end
L_DOP(DIP, DOP, POP, p) = -remineralization(DOP, p) + dissolution(POP, p)

function NL_POP(DIP, DOP, POP, p)
    @unpack σ = p
    return (1 - σ) * uptake(DIP, p)
end
L_POP(DIP, DOP, POP, p) = -dissolution(POP, p)

NL_all = (NL_DIP, NL_DOP, NL_POP)
L_all = (L_DIP, L_DOP, L_POP)

F2, L, NL, ∇ₓF2, ∇ₓL, ∇ₓNL, T = split_state_function_and_Jacobian(T_all, L_all, NL_all, nb)

#===========================================
Tests
===========================================#

@testset "Biogeochemical functions" begin
    @testset "Particle Flux Divergence is conservative" begin
        @test norm(v) / norm(T_POP(p)' * v) > ustrip(upreferred(1u"Myr"))
    end
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n) .* exp.(cos.(collect(1:n)) / 10)
    @testset "State function" begin
        F₀ = fun(x, p)
        v_all = repeat(v, nt)
        @test F₀ isa Vector{Float64}
        @test size(F₀) == (n,)
        @test norm(v_all .* x) / norm(v_all' * F₀) > ustrip(upreferred(1u"Myr"))
    end
    @testset "Jacobian of the state function" begin
        ∇ₓF₀ = fun(Val{:jac}, x, p)
        @test ∇ₓF₀ isa SparseMatrixCSC
        @test size(∇ₓF₀) == (n, n)
    end
end


@testset "In-place functions" begin
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n) .* exp.(cos.(collect(1:n)) / 10)
    dx = similar(x)
    @testset "in-place F ≈ out-of-place F" begin
        @test fun(x, p) ≈ fun!(dx, x, p) rtol = 1e-10
        @test fun(Val{:jac}, x, p) ≈ fun!(Val{:jac}, x, p) rtol = 1e-10
    end
end

@testset "Linear/nonlinear split functions" begin
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n) .* exp.(cos.(collect(1:n)) / 10)
    @testset "split F ≈ F" begin
        @test F2(x, p) ≈ fun(x, p) rtol = 1e-10
    end
    @testset "split F ≈ L + NL" begin
        @test F2(x, p) ≈ L(x, p) + NL(x, p) - T(p) * x rtol = 1e-10
    end
    @testset "split ∇ₓF ≈ ∇ₓF" begin
        @test fun(Val{:jac}, x, p) ≈ ∇ₓF2(x, p) rtol = 1e-10
    end
    @testset "split ∇ₓF ≈ T + ∇ₓNL" begin
        @test ∇ₓF2(x, p) ≈ ∇ₓL(p) + ∇ₓNL(x, p) - T(p) rtol = 1e-10
    end
end
