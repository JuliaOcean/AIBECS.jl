#---------------------------------------------------------
# # [P-cycle: 3 pools (DIP+DOP+POP)](@id P-cycle-3pool)
#---------------------------------------------------------

# This tutorial expands the [previous 2-pool DIP+POP P-cycle](@ref P-cycle-2pool)
# to also track dissolved organic phosphorus (DOP). Phosphorus is conserved
# exactly across the three pools: every mole removed by any uptake or
# dissolution/remineralization is allocated to one of DIP, DOP, or POP. A
# single shared fraction $\sigma$ partitions both uptake pathways (from DIP
# and from DOP) between DOP and POP, with the only basin-scale source/sink
# being the slow geological restoring on DIP.
#
# ```mermaid
# flowchart LR
#     geo([DIP_geo]) -- "1 / τ_geo" --> DIP
#     DIP -- "σ · U_DIP" --> DOP
#     DIP -- "(1 - σ) · U_DIP" --> POP
#     DOP -- "(1 - σ) · U_DOP" --> POP
#     POP -- "D_POP" --> DOP
#     DOP -- "R_DOP" --> DIP
# ```
#
# (The σ · U_DOP self-loop on DOP nets to zero and is omitted.)
#
# We use the OCCA circulation here (smaller than OCIM, lighter for CI docs).

using AIBECS
using JLD2                                        # required by `OCCA.load`
using Unitful
using Unitful: m, d, yr, Myr, mol, mmol, μmol, NoUnits
import AIBECS: @units, units
import AIBECS: @initial_value, initial_value

# --- Circulation -----------------------------------------------------------

grd, T_OCCA = OCCA.load()
T_DIP(p) = T_OCCA
T_DOP(p) = T_OCCA
T_POP(p) = transportoperator(grd, z -> w(z, p))

# --- Process rates ---------------------------------------------------------

function w(z, p)
    @unpack w₀, w′ = p
    return @. w₀ + w′ * z
end

z = depthvec(grd)

function U_DIP(DIP, p)
    @unpack τ_DIP, k_DIP, z₀ = p
    return @. DIP / τ_DIP * DIP / (DIP + k_DIP) * (z ≤ z₀) * (DIP ≥ 0)
end

function U_DOP(DOP, p)
    @unpack τ_DOPup, k_DOP, z₀ = p
    return @. DOP / τ_DOPup * DOP / (DOP + k_DOP) * (z ≤ z₀) * (DOP ≥ 0)
end

function R_DOP(DOP, p)
    @unpack τ_DOPrem = p
    return DOP / τ_DOPrem
end

function D_POP(POP, p)
    @unpack τ_POPdiss = p
    return POP / τ_POPdiss
end

# --- Source/sink (P-conserving) -------------------------------------------

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

# --- Parameters ------------------------------------------------------------

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

# --- Build and solve -------------------------------------------------------

nb = sum(iswet(grd))
F = AIBECSFunction((T_DIP, T_DOP, T_POP), (G_DIP, G_DOP, G_POP), nb)

@unpack DIP_geo = p
x0 = DIP_geo * ones(3nb)
prob = SteadyStateProblem(F, x0, p)
sol = solve(prob, CTKAlg(), preprint = " ").u

DIP, DOP, POP = state_to_tracers(sol, grd)

# --- Plot: single figure, 3 zonal-average panels --------------------------

using CairoMakie
using MakieExtra
using Unitful: ustrip
CairoMakie.activate!(type = "png")

# Disambiguate Makie vs. MakieExtra `get_ticks` for explicit (positions, labels)
# under a SymLog/AsinhScale colorbar — both packages provide a method that
# matches that combo, so we pin it: an explicit (positions, labels) tuple is
# always used verbatim, regardless of scale.
Makie.get_ticks(t::Tuple{Any, Any},
                ::Union{MakieExtra.SymLog, MakieExtra.AsinhScale},
                ::Makie.Automatic, vmin, vmax) = t

lat = ustrip.(grd.lat)
depth = ustrip.(grd.depth)
latticks = (-90:30:90, ["90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"])

zm_DIP = AIBECS.zonalmean(DIP * u"mol/m^3", grd, 1)
zm_DOP = AIBECS.zonalmean(DOP * u"mol/m^3", grd, 1)
zm_POP = AIBECS.zonalmean(POP * u"mol/m^3", grd, 1)

maxdepth = ceil(maximum(depth) / 1000) * 1000

# Snap hi/n to a 1/2/5 × 10^k "nice" step and return level boundaries starting at 0.
function nicelevels(data; n = 10)
    hi = maximum(filter(isfinite, vec(data)))
    raw = hi / n
    mag = 10.0 ^ floor(log10(raw))
    nice = (raw/mag < 1.5 ? 1.0 : raw/mag < 3.5 ? 2.0 : raw/mag < 7.5 ? 5.0 : 10.0) * mag
    return collect(0:nice:(ceil(hi / nice) * nice))
end

# Log-ish ladder: 0 plus multiples of {1, 2, 5} × 10^k spanning ~`decades`
# decades below max(data).
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

data_DIP = ustrip.(u"μmol/m^3", zm_DIP)
data_DOP = ustrip.(u"nM", zm_DOP)
data_POP = ustrip.(u"nM", zm_POP)

levels_DIP = nicelevels(data_DIP; n = 10)
levels_DOP = logishlevels(data_DOP)
levels_POP = logishlevels(data_POP)

panels = [
    (label = "A", name = "DIP", data = data_DIP,
     base_cmap = :viridis, levels = levels_DIP, scale = identity,
     cblabel = rich("DIP (µmol m", superscript("−3"), ")")),
    (label = "B", name = "DOP", data = data_DOP,
     base_cmap = :magma, levels = levels_DOP,
     scale = SymLog(symlog_thresh(levels_DOP)),
     cblabel = "DOP (nM)"),
    (label = "C", name = "POP", data = data_POP,
     base_cmap = :cividis, levels = levels_POP,
     scale = SymLog(symlog_thresh(levels_POP)),
     cblabel = "POP (nM)"),
]

fig = Figure(size = (820, 1000))
for (i, panel) in enumerate(panels)
    ax = Axis(fig[i, 1];
        backgroundcolor = :lightgray,
        xticks = latticks,
        ylabel = "Depth (m)",
        xlabel = i == 3 ? "Latitude" : "",
        yreversed = true,
        limits = (nothing, nothing, 0, maxdepth),
        xautolimitmargin = (0.0, 0.0),
        yautolimitmargin = (0.0, 0.0))
    i < 3 && hidexdecorations!(ax;
        ticklabels = true, label = true, ticks = false, grid = false)
    cmap = cgrad(panel.base_cmap, length(panel.levels) - 1; categorical = true)
    cf = CairoMakie.contourf!(ax, lat, depth, panel.data;
        levels = panel.levels, colormap = cmap, nan_color = :lightgray,
        extendlow = cmap[1], extendhigh = cmap[end],
        colorscale = panel.scale)
    if panel.scale === identity
        cbticks = panel.levels[1:2:end]
        Colorbar(fig[i, 2], cf;
            label = panel.cblabel,
            ticks = (cbticks, ticklabel.(cbticks)))
    else
        Colorbar(fig[i, 2], cf;
            label = panel.cblabel,
            ticks = BaseMulTicks([1]),
            minorticks = BaseMulTicks(2:9),
            minorticksvisible = true)
    end
    text!(ax, 0, 0; text = "$(panel.label). $(panel.name) zonal mean",
        space = :relative,
        align = (:left, :bottom), offset = (6, 6),
        font = :bold, color = :black, fontsize = 19)
end

fig
