using AIBECS
import AIBECS: @units, units
import AIBECS: @initial_value, initial_value
using Unitful
using Unitful: m, d, s, yr, Myr, mol, mmol, μmol, μM
using UnPack

# ---------------------------------------------------------------------------
# Tracer setups used by the benchmark harness. Each setup mirrors its
# corresponding tutorial in docs/lit/tutorials/ as closely as possible —
# same equations, same parameter values, same initial guess. The only
# difference is that the benchmark builders take (grd, T_circ) so they can
# be applied across the circulation tier matrix; tutorials hard-code a
# specific circulation.
# ---------------------------------------------------------------------------

# === Ideal age (mirrors 1_ideal_age.jl) ====================================

struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end

function build_idealage_problem(grd, T_circ)
    z = depthvec(grd)
    G(x, p) = (@unpack τ, z₀ = p; @. 1 - x / τ * (z ≤ z₀))
    F  = AIBECSFunction(T_circ, G)
    p  = IdealAgeParameters(1.0, 30.0)
    nb = sum(iswet(grd))
    return SteadyStateProblem(F, zeros(nb), p)
end

# === Radiocarbon (mirrors 2_radiocarbon.jl) ================================

@units struct RadiocarbonParameters{U} <: AbstractParameters{U}
    λ::U    | u"m/yr"
    h::U    | u"m"
    τ::U    | u"yr"
    Ratm::U | u"M"
end

function build_radiocarbon_problem(grd, T_circ)
    z = depthvec(grd)
    G(R, p) = (@unpack λ, h, Ratm, τ = p; @. λ / h * (Ratm - R) * (z ≤ h) - R / τ)
    F = AIBECSFunction(T_circ, G)
    p = RadiocarbonParameters(λ    = 50u"m"/10u"yr",
                              h    = grd.δdepth[1],
                              τ    = 5730u"yr"/log(2),
                              Ratm = 42.0u"nM")
    nb = sum(iswet(grd))
    return SteadyStateProblem(F, zeros(nb), p)
end

# === Coupled PO4-POP (mirrors 3_Pmodel.jl) =================================

@initial_value @units struct PmodelParameters{U} <: AbstractParameters{U}
    w₀::U       |  0.64 | m/d
    w′::U       |  0.13 | m/d/m
    τ_DIP::U    | 230.0 | d
    k::U        |  6.62 | μmol/m^3
    z₀::U       |  80.0 | m
    τ_POP::U    |   5.0 | d
    τ_geo::U    |   1.0 | Myr
    DIP_geo::U  |  2.12 | mmol/m^3
end

function build_pmodel_problem(grd, T_circ)
    nb = sum(iswet(grd))
    z  = depthvec(grd)

    T_DIP(p) = T_circ
    function w(z, p)
        @unpack w₀, w′ = p
        return @. w₀ + w′ * z
    end
    T_POP(p) = transportoperator(grd, z -> w(z, p))

    function U(x, p)
        @unpack τ_DIP, k, z₀ = p
        return @. x / τ_DIP * x / (x + k) * (z ≤ z₀) * (x ≥ 0)
    end
    function R(x, p)
        @unpack τ_POP = p
        return x / τ_POP
    end
    function G_DIP(DIP, POP, p)
        @unpack DIP_geo, τ_geo = p
        return @. -$U(DIP, p) + $R(POP, p) + (DIP_geo - DIP) / τ_geo
    end
    function G_POP(DIP, POP, p)
        @unpack τ_geo = p
        return @. $U(DIP, p) - $R(POP, p) - POP / τ_geo
    end

    F = AIBECSFunction((T_DIP, T_POP), (G_DIP, G_POP), nb)
    p = PmodelParameters()
    @unpack DIP_geo = p
    x = DIP_geo * ones(2nb)
    return SteadyStateProblem(F, x, p)
end

const TRACER_SETUPS = (
    idealage    = (label = "idealage",    build = build_idealage_problem),
    radiocarbon = (label = "radiocarbon", build = build_radiocarbon_problem),
    po4pop      = (label = "po4pop",      build = build_pmodel_problem),
)
