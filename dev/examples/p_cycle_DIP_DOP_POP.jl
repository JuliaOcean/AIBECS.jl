
#---------------------------------------------------------
# # A Phosphorus Cycling Model
#---------------------------------------------------------

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`p_cycle_DIP_DOP_POP.ipynb`](@__NBVIEWER_ROOT_URL__examples/generated/p_cycle_DIP_DOP_POP.ipynb)

#---------------------------------------------------------
# ## Tracer equations
#---------------------------------------------------------

# We consider a simple model for the cycling of phosphorus with 3 state variables consisting of phosphate (PO₄) AKA dissolved inorganic phosphorus (DIP), dissolved organic phosphorus (DOP), and particulate organic phosphorus (POP).
#
# The dissolved phases are transported by advection and diffusion whereas the particulate phase sinks rapidly down the water column without any appreciable transport by the circulation.
#
# The governing equations are:
#
# $$\frac{\partial}{\partial t} DIP + \nabla \cdot \left[\boldsymbol{u} + \mathbf{K}\cdot\nabla \right] DIP = -\gamma(DIP) + \kappa_\mathsf{D} \, DOP,$$
#
# $$\frac{\partial}{\partial t} DOP + \nabla \cdot \left[\boldsymbol{u} + \mathbf{K}\cdot \nabla \right] DOP = \sigma \, \gamma(DIP) + \kappa_\mathsf{P} \, POP - \kappa_\mathsf{D} \, DOP,$$
#
# and
#
# $$\frac{\partial}{\partial t} POP + \frac{\partial}{\partial z} \left[w_\mathsf{P} \, POP\right] = (1-\sigma) \, \gamma(DIP) - \kappa_\mathsf{P} \, POP,$$
#
# where $\boldsymbol{u}$ is the fluid velocity and $\mathbf{K}$ is the eddy-diffusion tensor.
# Thus, $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right]$ is a differential operator that represents the transport by the ocean circulation.
# The function $\gamma(DIP)$ represents the biological uptake of DIP by phytoplankton.

# Oxygen participates to this cycle too and satisfies its own tracer equation
#
# $$\frac{\partial}{\partial t} DO_2 + \nabla \cdot \left[\boldsymbol{u} + \mathbf{K}\cdot\nabla \right] DO_2 = -r_{\mathsf{O}_2:\mathsf{P}} \, \kappa_\mathsf{D} \, DOP + \boldsymbol{\Lambda}(DO_2 - [O_2]_{\mathsf{sat}})$$
#
# where $\Lambda$ is the air-sea gas exchange operator.

# These tracer equations depend on a number of scalars, that we list below
#
# | Symbol                        | Definition                                                                  |
# |:------------------------------|:----------------------------------------------------------------------------|
# | $w_\mathsf{P}$                | depth dependent particle sinking speed                                      |
# | $\sigma$                      | fraction of the organic matter production allocated to the dissolved phase  |
# | $\kappa_\mathsf{D}$           | respiration rate for dissolved organic matter (DOP → DIP)                   |
# | $\kappa_\mathsf{P}$           | dissolution rate for particulate organic matter (POP → DOP)                 |
# | $r_{\mathsf{O}_2:\mathsf{P}}$ | number of moles of O₂ needed to respire 1 mole of DOP                       |

using AIBECS

# Load the circulation and grid

const wet3d, grd, T_Circulation = OCIM1.load()

# Define useful constants and arrays

const iwet = indices_of_wet_boxes(wet3d)
const nb = number_of_wet_boxes(wet3d)
const v = vector_of_volumes(wet3d, grd)
const z = vector_of_depths(wet3d, grd)
const ztop = vector_of_top_depths(wet3d, grd)

# And matrices

const DIV = buildDIV(wet3d, iwet, grd)
const Iabove = buildIabove(wet3d, iwet) ;

# ### Transport matrices

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

# Because AIBECS will solve for the steady state solution directly without time-stepping the goverining equations to equilibrium, we don't have any opportunity to specify any intial conditions.
# Initial conditions are how the total amount of conserved elements get specified in most global biogeochemical modelels.
# Thus to specify the total inventory of P in AIBECS we add a very weak resporing term to the DIP equation.
# The time-scale for this restoring term is chosen to be very long compared to the timescale with which the ocean circulation homogenizes a tracer.
# Because of this long timescale we call it the geological restoring term, but geochemists who work on geological processes don't like that name!
# In any event the long timescale allows us to prescribe the total inventory of P without having any appreciable impact on the 3d distribution of P.

# ### Sources minus sinks
#
# ##### Geological Restoring

function geores(x, p)
    τg, xgeo = p.τg, p.xgeo
    return (xgeo .- x) / τg
end

# ##### Uptake of phosphate (DIP)

relu(x) = (x .≥ 0) .* x
function uptake(DIP, p)
    τu, ku, z₀ = p.τu, p.ku, p.z₀
    DIP⁺ = relu(DIP)
    return 1/τu * DIP⁺.^2 ./ (DIP⁺ .+ ku) .* (z .≤ z₀)
end

# ##### Remineralization DOP into DIP

function remineralization(DOP, p)
    κDOP = p.κDOP
    return κDOP * DOP
end

# ##### Dissolution of POP into DOP

function dissolution(POP, p)
    κPOP = p.κPOP
    return κPOP * POP
end

# ##### Air-sea gas exchange

dz1 = grd["dzt"][1]               # thickness of the top layer
z = vec(grd["ZT3d"])[iwet]        # depth of the gridbox centers
using WorldOceanAtlasTools
WOA = WorldOceanAtlasTools
μDO2 , σ²DO2 = WOA.fit_to_grid(grd,2018,"O2sat","annual","1°","an")

# This below was commented out?

function airsea(DO2, p)
    κDO2 = p.κDO2
    return κDO2 * (z .< 20) .* (1.0 .- DO2) / dz1
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


# ### Parameters
#
# Build the parameters type and p₀

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

# ### Generate state function and Jacobian

nt = length(T_all)    # number of tracers
n = nt * nb           # total dimension of the state vector
p = Pcycle_Parameters() # parameters
x = p.xgeo * ones(n) # initial iterate
# F, ∇ₓF = state_function_and_Jacobian(T_all, sms_all, nb)

# and solve

# prob = SteadyStateProblem(F, ∇ₓF, x, p)
# s = solve(prob, CTKAlg());

# unpack state

# function unpack_P_state(s,mask)
#         DIP = NaN*mask; DOP = NaN*mask; POP = NaN*mask
#         iwet = findall(x-> x==1, vec(mask))
#         nwet = length(iwet)
#         idip = 1:nwet;       idop = idip.+nwet;   ipop = idop.+nwet
#         DIP[iwet] = s[idip]; DOP[iwet] = s[idop]; POP[iwet] = s[ipop]
#     return DIP, DOP, POP
# end
# DIP, DOP, POP = unpack_P_state(s,wet3d);

# We will plot the concentration of DIP at a given depth horizon

depth = vec(grd["zt"])
iz = findfirst(depth .> 200)
iz, depth[iz]

#-

# dip = DIP[:,:,iz] * ustrip(1.0u"mol/m^3"|>u"mmol/m^3");
# lat, lon = vec(grd["yt"]), vec(grd["xt"]);

# and plot

# ENV["MPLBACKEND"]="qt5agg"
# using PyPlot, PyCall
# using Conda; Conda.add("Cartopy")
# clf()
# ccrs = pyimport("cartopy.crs")
# ax = subplot(projection=ccrs.Robinson(central_longitude=-155.0))
# ax.coastlines()
# # making it cyclic for Cartopy
# lon_cyc = [lon; 360+lon[1]]
# dip_cyc = hcat(dip, dip[:,1])
# # And plot
# p = contourf(lon_cyc, lat, dip_cyc, levels=0:0.2:3.6, transform=ccrs.PlateCarree(), zorder=-1)
# colorbar(p, orientation="horizontal");
# gcf()
