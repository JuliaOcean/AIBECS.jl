#---------------------------------------------------------
# # [P-cycle: 3 pools (DIP+DOP+POP)](@id P-cycle-3pool)
#---------------------------------------------------------

# This tutorial expands the [previous 2-pool DIP+POP P-cycle](@ref P-cycle-2pool)
# to also track dissolved organic phosphorus (DOP). Phosphorus is conserved
# exactly across the three pools: every mole removed by any uptake or
# dissolution/remineralization is allocated to one of DIP, DOP, or POP. Both
# DIP and DOP are taken up by phytoplankton, and a shared fraction $\sigma$
# of each uptake is routed to DOP with the remaining $(1 - \sigma)$ routed
# to POP. The only basin-scale source/sink is a slow geological restoring
# on DIP.

# We use the OCIM0 circulation.

using AIBECS
using JLD2                                        # required by `OCIM0.load`
using Unitful
using Unitful: m, d, yr, Myr, mol, mmol, μmol, NoUnits
import AIBECS: @units, units
import AIBECS: @initial_value, initial_value

# The dissolved phases (DIP, DOP) ride the OCIM0 circulation; the particulate
# phase (POP) sinks with depth-dependent speed $w(z) = w_0 + w' z$.

grd, T_OCIM0 = OCIM0.load()
T_DIP(p) = T_OCIM0
T_DOP(p) = T_OCIM0
T_POP(p) = transportoperator(grd, z -> w(z, p))

function w(z, p)
    @unpack w₀, w′ = p
    return @. w₀ + w′ * z
end

z = depthvec(grd)

# Uptake by phytoplankton (DIP and DOP), saturating in concentration and
# restricted to the euphotic layer $z \le z_0$:

function U_DIP(DIP, p)
    @unpack τ_DIP, k_DIP, z₀ = p
    return @. DIP / τ_DIP * DIP / (DIP + k_DIP) * (z ≤ z₀) * (DIP ≥ 0)
end

function U_DOP(DOP, p)
    @unpack τ_DOPup, k_DOP, z₀ = p
    return @. DOP / τ_DOPup * DOP / (DOP + k_DOP) * (z ≤ z₀) * (DOP ≥ 0)
end

# Remineralization (DOP → DIP) and dissolution (POP → DOP), linear:

function R_DOP(DOP, p)
    @unpack τ_DOPrem = p
    return DOP / τ_DOPrem
end

function D_POP(POP, p)
    @unpack τ_POPdiss = p
    return POP / τ_POPdiss
end

# The per-pool source/sink terms — note the shared $\sigma$ in the uptake
# partition, and the geological restoring on DIP:

function G_DIP(DIP, DOP, POP, p)
    @unpack DIP_geo, τ_geo = p
    return @. -$U_DIP(DIP, p) + $R_DOP(DOP, p) + (DIP_geo - DIP) / τ_geo
end

function G_DOP(DIP, DOP, POP, p)
    @unpack σ = p
    return @. σ * $U_DIP(DIP, p) - $R_DOP(DOP, p) - (1 - σ) * $U_DOP(DOP, p) + $D_POP(POP, p)
end

function G_POP(DIP, DOP, POP, p)
    @unpack σ = p
    return @. (1 - σ) * $U_DIP(DIP, p) + (1 - σ) * $U_DOP(DOP, p) - $D_POP(POP, p)
end

# Parameters:

@initial_value @units struct PmodelParameters{U} <: AbstractParameters{U}
    w₀::U        | 0.64   | m / d
    w′::U        | 0.13   | m / d / m
    τ_DIP::U     | 230.0  | d
    k_DIP::U     | 6.62   | μmol / m^3
    τ_DOPup::U   | 350.0  | d
    k_DOP::U     | 5.0    | μmol / m^3
    z₀::U        | 80.0   | m
    τ_DOPrem::U  | 0.25   | yr
    τ_POPdiss::U | 5.0    | d
    σ::U         | 1 / 3  | NoUnits
    τ_geo::U     | 1.0    | Myr
    DIP_geo::U   | 2.12   | mmol / m^3
end

p = PmodelParameters()

# Build and solve:

nb = sum(iswet(grd))
F = AIBECSFunction((T_DIP, T_DOP, T_POP), (G_DIP, G_DOP, G_POP), nb)
@unpack DIP_geo = p
x0 = DIP_geo * ones(3nb)
prob = SteadyStateProblem(F, x0, p)
sol = solve(prob, CTKAlg(), preprint = " ")

# The `sol` value above carries the retcode, stats, and the residual norm —
# use it to confirm convergence before unpacking the state.

DIP, DOP, POP = state_to_tracers(sol.u, grd)

# We plot zonal averages of the 3 tracers (rows) across the Atlantic,
# Pacific, and Indian basins (columns), with one colorbar per row on
# the right.

using CairoMakie
using MakieExtra
using OceanBasins
using Unitful: ustrip
CairoMakie.activate!(type = "png")

Makie.get_ticks(t::Tuple{Any, Any},
                ::Union{MakieExtra.SymLog, MakieExtra.AsinhScale},
                ::Makie.Automatic, vmin, vmax) = t

OCEANS = oceanpolygons()
masks = (
    Atlantic = isatlantic(latvec(grd), lonvec(grd), OCEANS),
    Pacific  = ispacific( latvec(grd), lonvec(grd), OCEANS),
    Indian   = isindian(  latvec(grd), lonvec(grd), OCEANS),
)
basins = (:Atlantic, :Pacific, :Indian)

lat = ustrip.(grd.lat)
depth = ustrip.(grd.depth)
latticks = (-90:30:90, ["90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"])
maxdepth = ceil(maximum(depth) / 1000) * 1000

function wetlatrange(zm, pad = 5)
    haswet = [any(!isnan, view(zm, i, :)) for i in eachindex(lat)]
    wetidx = findall(haswet)
    isempty(wetidx) && return (-90.0, 90.0)
    return (max(-90, lat[first(wetidx)] - pad),
            min( 90, lat[last(wetidx)]  + pad))
end

function nicelevels(data; n = 10)
    hi = maximum(filter(isfinite, vec(data)))
    raw = hi / n
    mag = 10.0 ^ floor(log10(raw))
    nice = (raw/mag < 1.5 ? 1.0 : raw/mag < 3.5 ? 2.0 : raw/mag < 7.5 ? 5.0 : 10.0) * mag
    return collect(0:nice:(ceil(hi / nice) * nice))
end

function logishlevels(data; decades = 3)
    hi = maximum(filter(isfinite, vec(data)))
    hi_pow = Int(ceil(log10(hi)))
    lo_pow = hi_pow - decades
    candidates = sort!(vcat([0.0],
        [a * 10.0^k for k in lo_pow:hi_pow for a in (1.0, 2.0, 5.0)]))
    keep = candidates[candidates .≤ hi]
    if !(hi in keep)
        above = findfirst(c -> c > hi, candidates)
        above !== nothing && push!(keep, candidates[above])
    end
    return keep
end

ticklabel(x) = isinteger(x) ? string(Int(x)) : string(round(x; sigdigits = 3))
symlog_thresh(levels) = (nz = filter(>(0), levels); isempty(nz) ? 1.0 : minimum(nz))

zm_data = Dict(
    (:DIP, b) => ustrip.(u"μM", AIBECS.zonalmean(DIP * u"mol/m^3", grd, masks[b])) for b in basins
)
merge!(zm_data, Dict(
    (:DOP, b) => ustrip.(u"nM", AIBECS.zonalmean(DOP * u"mol/m^3", grd, masks[b])) for b in basins
))
merge!(zm_data, Dict(
    (:POP, b) => ustrip.(u"nM", AIBECS.zonalmean(POP * u"mol/m^3", grd, masks[b])) for b in basins
))

allvals(t) = vcat((vec(zm_data[(t, b)]) for b in basins)...)
levels_DIP = nicelevels(allvals(:DIP); n = 10)
levels_DOP = logishlevels(allvals(:DOP))
levels_POP = logishlevels(allvals(:POP))

xlim_basin = Dict(b => wetlatrange(zm_data[(:DIP, b)]) for b in basins)

rows = [
    (name = :DIP, levels = levels_DIP, base_cmap = :viridis, scale = identity,
     cblabel = rich("DIP (µM)")),
    (name = :DOP, levels = levels_DOP, base_cmap = :magma,
     scale = SymLog(symlog_thresh(levels_DOP)), cblabel = "DOP (nM)"),
    (name = :POP, levels = levels_POP, base_cmap = :cividis,
     scale = SymLog(symlog_thresh(levels_POP)), cblabel = "POP (nM)"),
]

fig = Figure(size = (1400, 950))
for (i, row) in enumerate(rows)
    cmap = cgrad(row.base_cmap, length(row.levels) - 1; categorical = true)
    cf_last = nothing
    for (j, b) in enumerate(basins)
        ax = Axis(fig[i, j];
            backgroundcolor = :lightgray,
            xticks = latticks,
            yreversed = true,
            limits = (xlim_basin[b][1], xlim_basin[b][2], 0, maxdepth),
            xautolimitmargin = (0.0, 0.0),
            yautolimitmargin = (0.0, 0.0),
            xlabel = i == 3 ? "Latitude" : "",
            ylabel = j == 1 ? "Depth (m)" : "",
            title = i == 1 ? string(b) : "")
        i < 3 && hidexdecorations!(ax;
            ticklabels = true, label = true, ticks = false, grid = false)
        j > 1 && hideydecorations!(ax;
            ticklabels = true, label = true, ticks = false, grid = false)
        cf_last = CairoMakie.contourf!(ax, lat, depth, zm_data[(row.name, b)];
            levels = row.levels, colormap = cmap, nan_color = :lightgray,
            extendlow = cmap[1], extendhigh = cmap[end],
            colorscale = row.scale)
    end
    if row.scale === identity
        cbticks = row.levels[1:2:end]
        Colorbar(fig[i, length(basins) + 1], cf_last;
            label = row.cblabel,
            ticks = (cbticks, ticklabel.(cbticks)))
    else
        Colorbar(fig[i, length(basins) + 1], cf_last;
            label = row.cblabel,
            ticks = BaseMulTicks([1]),
            minorticks = BaseMulTicks(2:9),
            minorticksvisible = true)
    end
end

for (j, b) in enumerate(basins)
    colsize!(fig.layout, j, Auto(xlim_basin[b][2] - xlim_basin[b][1]))
end

fig
