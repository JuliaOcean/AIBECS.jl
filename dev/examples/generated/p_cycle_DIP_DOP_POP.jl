using AIBECS

const wet3d, grd, T_Circulation = OCIM1.load()

const iwet = indices_of_wet_boxes(wet3d)
const nb = number_of_wet_boxes(wet3d)
const v = vector_of_volumes(wet3d, grd)
const z = vector_of_depths(wet3d, grd)
const ztop = vector_of_top_depths(wet3d, grd)

const DIV = buildDIV(wet3d, iwet, grd)
const Iabove = buildIabove(wet3d, iwet) ;

T_DIP(p) = T_Circulation
T_DOP(p) = T_Circulation
T_DO2(p) = T_Circulation
const S₀ = buildPFD(ones(nb), DIV, Iabove)
const S′ = buildPFD(ztop, DIV, Iabove)
function T_POP(p)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S₀ + w′ * S′
end
T_all = (T_DIP, T_DOP, T_POP, T_DO2)

function geores(x, p)
    τg, xgeo = p.τg, p.xgeo
    return (xgeo .- x) / τg
end

relu(x) = (x .≥ 0) .* x
function uptake(DIP, p)
    τu, ku, z₀ = p.τu, p.ku, p.z₀
    DIP⁺ = relu(DIP)
    return 1/τu * DIP⁺.^2 ./ (DIP⁺ .+ ku) .* (z .≤ z₀)
end

function remineralization(DOP, p)
    κDOP = p.κDOP
    return κDOP * DOP
end

function dissolution(POP, p)
    κPOP = p.κPOP
    return κPOP * POP
end

dz1 = grd["dzt"][1]               # thickness of the top layer
z = vec(grd["ZT3d"])[iwet]        # depth of the gridbox centers
using WorldOceanAtlasTools
WOA = WorldOceanAtlasTools
μDO2 , σ²DO2 = WOA.fit_to_grid(grd,2018,"O2sat","annual","1°","an")

function airsea(DO2, p)
    κDO2 = p.κDO2
    return κDO2 * (z .< 20) .* (1.0 .- DO2) / dz1
end

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

t = empty_parameter_table()    # initialize table of parameters
add_parameter!(t, :xgeo, 2.17u"mmol/m^3",
    variance_obs = ustrip(upreferred(0.1 * 2.17u"mmol/m^3"))^2,
    description = "Geological mean P concentration",
    LaTeX = "\\state^\\mathrm{geo}")
add_parameter!(t, :τg, 1.0u"Myr",
    description = "Geological restoring timescale",
    LaTeX = "\\tau_\\mathrm{geo}")
add_parameter!(t, :ku, 10.0u"μmol/m^3",
    optimizable = true,
    description = "Half-saturation constant (Michaelis-Menten)",
    LaTeX = "k_\\vec{u}")
add_parameter!(t, :z₀, 80.0u"m",
    description = "Depth of the euphotic layer base",
    LaTeX = "z_0")
add_parameter!(t, :w₀, 1.0u"m/d",
    optimizable = true,
    description = "Sinking velocity at surface",
    LaTeX = "w_0")
add_parameter!(t, :w′, 1/4.4625u"d",
    optimizable = true,
    description = "Vertical gradient of sinking velocity",
    LaTeX = "w'")
add_parameter!(t, :κDOP, 1/0.25u"yr",
    optimizable = true,
    description = "Remineralization rate constant (DOP to DIP)",
    LaTeX = "\\kappa")
add_parameter!(t, :κPOP, 1/5.25u"d",
    optimizable = true,
    description = "Dissolution rate constant (POP to DOP)",
    LaTeX = "\\kappa")
add_parameter!(t, :σ, 0.3u"1",
    description = "Fraction of quick local uptake recycling",
    LaTeX = "\\sigma")
add_parameter!(t, :τu, 30.0u"d",
    optimizable = true,
    description = "Maximum uptake rate timescale",
    LaTeX = "\\tau_\\vec{u}")
initialize_Parameters_type(t, "Pcycle_Parameters")   # Generate the parameter type

nt = length(T_all)    # number of tracers
n = nt * nb           # total dimension of the state vector
p = Pcycle_Parameters() # parameters
x = p.xgeo * ones(n) # initial iterate

depth = vec(grd["zt"])
iz = findfirst(depth .> 200)
iz, depth[iz]

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

