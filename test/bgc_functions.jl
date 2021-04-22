#===========================================
Transport matrices
===========================================#
T_D = T_Circulation
function w(p)
    @unpack w₀, w′ = p
    z -> w₀ + w′ * z
end
T_POP = iiptransportoperator(grd; w)
update_coefficients!(T_POP, nothing, p, nothing)
T_POP_old(p) = transportoperator(grd, w(p))
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
Gs = (sms_DIP, sms_DOP, sms_POP) # bundles all the source-sink functions in a tuple

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
G!s = (G_DIP!, G_DOP!, G_POP!)


#===========================================
AIBECS F and ∇ₓF
===========================================#
fun = AIBECSFunction(T_all, Gs)
F, ∇ₓF = F_and_∇ₓF(T_all, Gs)

fun1 = AIBECSFunction(T_all, Gs)
fun2 = AIBECSFunction(T_all, G!s)
F1, ∇ₓF1 = F_and_∇ₓF(T_all, Gs)
F2, ∇ₓF2 = F_and_∇ₓF(T_all, G!s)

#===========================================
AIBECS split operators for CNLF
===========================================#
# TODO Remove deprecated split functionality

NL_DIP!(dx, DIP, DOP, POP, p) = @. dx .= -$uptake(DIP, p) + $geores(0*DIP, p)
#L_DIP(DIP, DOP, POP, p) = remineralization(DOP, p) + geores(DIP, p) - geores(0*DIP, p)

function NL_DOP!(dx, DIP, DOP, POP, p)
    @unpack σ = p
    @. dx = σ * $uptake(DIP, p)
end
#L_DOP(DIP, DOP, POP, p) = -remineralization(DOP, p) + dissolution(POP, p)

function NL_POP!(dx, DIP, DOP, POP, p)
    @unpack σ = p
    @. dx = (1 - σ) * $uptake(DIP, p)
end
#L_POP(DIP, DOP, POP, p) = -dissolution(POP, p)

NL_all! = (NL_DIP!, NL_DOP!, NL_POP!)
#_all = (L_DIP, L_DOP, L_POP)

#F3, L, NL, ∇ₓF3, ∇ₓL, ∇ₓNL, T = split_state_function_and_Jacobian(T_all, L_all, NL_all, nb)



#==========================================================
AIBECS split operators with new in-place operator framework
==========================================================#
nt = 3
@enum Tracers dip=1 dop pop
function update_remin_op(oldval,u,p,t)
    @unpack τDOP = p
    1 / τDOP
end
remin_op = SubBlockOperator(DiffEqScalar(1.0, update_func=update_remin_op), Int(dip), Int(dop), Int(dop))
function update_diss_op(oldval,u,p,t)
    @unpack τPOP = p
    1 / τPOP
end
diss_op = SubBlockOperator(DiffEqScalar(1.0, update_func=update_diss_op), Int(dop), Int(pop), Int(pop))
function geores_const(p)
    @unpack τgeo, xgeo = p
    xgeo / τgeo
end
function update_geores_op(oldval,u,p,t)
    @unpack τgeo = p
    1 / τgeo
end
geores_op = SubBlockOperator(DiffEqScalar(1.0, update_func=update_geores_op), nothing, Int(dip), Int(dip))
L_all = (remin_op, diss_op, geores_op)


split_fun = split_AIBECSFunction(T_all, L_all, NL_all!, nb)


#===========================================
Tests
===========================================#

@testset "Biogeochemical functions" begin
    @test T_POP_old(p) ≈ T_POP.A
    @testset "Particle Flux Divergence is conservative" begin
        @test norm(v) / norm(T_POP.A' * v) > ustrip(upreferred(1u"Myr"))
    end
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n) .* exp.(cos.(collect(1:n))/10)
    @testset "State function" begin
        F₀ = F1(x, p)
        v_all = repeat(v, nt)
        @test F₀ isa Vector{Float64}
        @test size(F₀) == (n,)
        # Not sure about this test below, which sort-of tests if F is conservative,
        # although it should not be because of the geological restoring.
        @test norm(v_all .* x) / norm(v_all' * F₀) > ustrip(upreferred(1u"Myr"))
    end
    @testset "Jacobian of the state function" begin
        ∇ₓF₀ = ∇ₓF1(x, p)
        @test ∇ₓF₀ isa SparseMatrixCSC
        @test size(∇ₓF₀) == (n,n)
    end
end


@testset "In-place functions" begin
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n) .* exp.(cos.(collect(1:n))/10)
    dx = similar(x)
    @testset "F! ≈ F" begin
        F₀ = F2(dx, x, p)
        @test F₀ ≈ F1(x,p) rtol=1e-10
    end
end

# TODO: Rewrite tests for new split API
#@testset "Linear/nonlinear split functions" begin
#    nt = length(T_all)
#    n = nt * nb
#    @unpack xgeo = p
#    x = xgeo * ones(n) .* exp.(cos.(collect(1:n))/10)
#    @testset "split F ≈ F" begin
#        @test F3(x,p) ≈ F1(x,p) rtol=1e-10
#    end
#    @testset "split F ≈ L + NL" begin
#        @test F3(x,p) ≈ L(x,p) + NL(x,p) - T(p) * x rtol=1e-10
#    end
#    @testset "split ∇ₓF ≈ ∇ₓF" begin
#        @test ∇ₓF1(x,p) ≈ ∇ₓF3(x,p) rtol=1e-10
#    end
#    @testset "split ∇ₓF ≈ T + ∇ₓNL" begin
#        @test ∇ₓF3(x,p) ≈ ∇ₓL(p) + ∇ₓNL(x,p) - T(p) rtol=1e-10
#    end
#end
