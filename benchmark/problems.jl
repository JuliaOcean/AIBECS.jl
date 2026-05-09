# Benchmark problem builders, shared between the benchmark harness
# (benchmark/solvers.jl) and the test suite (test/sparse_jacobian.jl).
#
# Each builder returns a `SteadyStateProblem` that mirrors a tutorial
# (docs/lit/tutorials/{1_ideal_age,2_radiocarbon,3_Pmodel}.jl).

using AIBECS
using UnPack
using Unitful
using Unitful: m, d, yr, Myr, mmol, μmol
using SciMLBase
import AIBECS: @units, units
import AIBECS: @initial_value, initial_value

# Ideal age — mirrors docs/lit/tutorials/1_ideal_age.jl
struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end
function build_idealage_problem(grd, T_circ)
    z = depthvec(grd)
    G(x, p) = (@unpack τ, z₀ = p; @. 1 - x / τ * (z ≤ z₀))
    F = AIBECSFunction(T_circ, G)
    p = IdealAgeParameters(1.0, z[1])
    nb = sum(iswet(grd))
    return SteadyStateProblem(F, zeros(nb), p)
end

# Radiocarbon — mirrors docs/lit/tutorials/2_radiocarbon.jl
@units struct RadiocarbonParameters{U} <: AbstractParameters{U}
    λ::U | u"m/yr"
    h::U | u"m"
    τ::U | u"yr"
    Ratm::U | u"M"
end
function build_radiocarbon_problem(grd, T_circ)
    z = depthvec(grd)
    G(R, p) = (@unpack λ, h, Ratm, τ = p; @. λ / h * (Ratm - R) * (z ≤ h) - R / τ)
    F = AIBECSFunction(T_circ, G)
    p = RadiocarbonParameters(
        λ = 50u"m" / 10u"yr",
        h = grd.δdepth[1],
        τ = 5730u"yr" / log(2),
        Ratm = 42.0u"nM",
    )
    nb = sum(iswet(grd))
    return SteadyStateProblem(F, zeros(nb), p)
end

# Coupled PO4–POP — mirrors docs/lit/tutorials/3_Pmodel.jl
@initial_value @units struct PmodelParameters{U} <: AbstractParameters{U}
    w₀::U | 0.64 | m / d
    w′::U | 0.13 | m / d / m
    τ_DIP::U | 230.0 | d
    k::U | 6.62 | μmol / m^3
    z₀::U | 80.0 | m
    τ_POP::U | 5.0 | d
    τ_geo::U | 1.0 | Myr
    DIP_geo::U | 2.12 | mmol / m^3
end
function build_pmodel_problem(grd, T_circ)
    nb = sum(iswet(grd))
    z = depthvec(grd)
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
    @unpack DIP_geo, τ_DIP, τ_POP = p
    # Initial guess at each tracer's physical scale: DIP at DIP_geo, POP at
    # the surface-mass-balance estimate POP_scale ≈ DIP_geo · τ_POP/τ_DIP
    # (matches `pmodel_scales`).
    POP_scale = DIP_geo * (τ_POP / τ_DIP)
    u0 = vcat(fill(DIP_geo, nb), fill(POP_scale, nb))
    return SteadyStateProblem(F, u0, p)
end

# Physics-based scale recipe for the coupled DIP/POP system. AIBECS convention
# (`src/Parameters.jl:287`): `@unpack` on an `APar` returns
# `ustrip(upreferred(value * units))` — i.e. SI base units, dimensionless. So
# the values we read here are numerically what `f(u, p)` and `prob.f.jac`
# actually compute on.
#
# `scale_u_POP ≈ DIP_geo · τ_POP / τ_DIP` from the surface-box mass balance
# `U(DIP) = R(POP)` once `DIP_geo ≫ k`; `scale_F` is the matching local rate
# `scale_u / τ_characteristic`.
function pmodel_scales(prob)
    @unpack DIP_geo, τ_DIP, τ_POP = prob.p     # SI: mol/m³, s, s
    nb = length(prob.u0) ÷ 2
    DIP_scale   = DIP_geo
    POP_scale   = DIP_scale * (τ_POP / τ_DIP)
    F_DIP_scale = DIP_scale / τ_DIP
    F_POP_scale = POP_scale / τ_POP            # = F_DIP_scale by construction
    return (
        vcat(fill(DIP_scale, nb),   fill(POP_scale, nb)),
        vcat(fill(F_DIP_scale, nb), fill(F_POP_scale, nb)),
    )
end
