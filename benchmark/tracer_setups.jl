using AIBECS
using AIBECS: @units, units
using Unitful
using UnPack

# ---------------------------------------------------------------------------
# Tracer setups used by the benchmark harness. Each setup returns a NamedTuple
# (label, build_problem) where build_problem(grd, T) constructs and returns a
# SteadyStateProblem ready to solve.
# ---------------------------------------------------------------------------

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

# Coupled DIP/POP P-cycle, mirroring docs/lit/tutorials/3_Pmodel.jl in shape
@units struct PmodelParameters{U} <: AbstractParameters{U}
    xgeo::U  | u"mmol/m^3"
    τgeo::U  | u"Myr"
    k::U     | u"μmol/m^3"
    w₀::U    | u"m/d"
    w′::U    | u"d^-1"
    τPOP::U  | u"d"
    τDIP::U  | u"d"
end

function build_pmodel_problem(grd, T_circ)
    nb     = sum(iswet(grd))
    iwet   = indices_of_wet_boxes(grd)
    ztop   = topdepthvec(grd)
    DIV    = DIVO(grd)
    Iabove = buildIabove(grd)
    S₀     = PFDO(ones(nb), DIV, Iabove)
    S′     = PFDO(ztop, DIV, Iabove)

    T_DIP(p) = T_circ
    function T_POP(p)
        @unpack w₀, w′ = p
        return w₀ * S₀ + w′ * S′
    end

    relu(x) = (x .≥ 0) .* x
    function uptake(DIP, p)
        @unpack τDIP, k = p
        DIP⁺ = relu(DIP)
        return @. 1 / τDIP * DIP⁺^2 / (DIP⁺ + k) * (ztop == 0)
    end
    function geores(x, p)
        @unpack τgeo, xgeo = p
        return @. (xgeo - x) / τgeo
    end

    sms_DIP(DIP, POP, p) = -uptake(DIP, p) .+ geores(DIP, p)
    sms_POP(DIP, POP, p) =  uptake(DIP, p) .- POP / 5.25u"d" |> ustrip

    F = AIBECSFunction((T_DIP, T_POP), (sms_DIP, sms_POP), nb)
    p = PmodelParameters(xgeo = 2.17u"mmol/m^3",
                         τgeo = 1.0u"Myr",
                         k    = 10.0u"μmol/m^3",
                         w₀   = 1.0u"m/d",
                         w′   = 0.22u"d^-1",
                         τPOP = 5.25u"d",
                         τDIP = 30.0u"d")
    x0 = ustrip(2.17e-3) * ones(2 * nb)   # uniform DIP initial guess
    return SteadyStateProblem(F, x0, p)
end

const TRACER_SETUPS = (
    idealage    = (label = "idealage",    build = build_idealage_problem),
    radiocarbon = (label = "radiocarbon", build = build_radiocarbon_problem),
    po4pop      = (label = "po4pop",      build = build_pmodel_problem),
)
