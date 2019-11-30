
#---------------------------------------------------------
# # A Phosphorus Cycling Model
#---------------------------------------------------------

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`p_cycle_DIP_DOP_POP.ipynb`](@__NBVIEWER_ROOT_URL__/examples/generated/p_cycle_DIP_DOP_POP.ipynb)

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
# $$\frac{\partial}{\partial t} DO_2 + \nabla \cdot \left[\boldsymbol{u} + \mathbf{K}\cdot\nabla \right] DO_2 = -r_{\mathsf{O}_2:\mathsf{P}} \, \kappa_\mathsf{D} \, DOP + \Lambda(DO_2 - [O_2]_{\mathsf{sat}})$$
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

grd, T_Circulation = OCIM0.load() ;

# Define useful constants and arrays

iwet = indices_of_wet_boxes(grd)
nb = length(iwet)
z = ustrip.(vector_of_depths(grd))
ztop = vector_of_top_depths(grd) ;

# And matrices

DIV = buildDIV(grd)
Iabove = buildIabove(grd) ;

# ### Transport matrices

T_DIP(p) = T_Circulation
T_DOP(p) = T_Circulation
T_DO2(p) = T_Circulation
S₀ = buildPFD(ones(nb), DIV, Iabove)
S′ = buildPFD(ztop, DIV, Iabove)
function T_POP(p)
    w₀, w′ = p.w₀, p.w′
    return w₀ * S₀ + w′ * S′
end
T_all = (T_DIP, T_DOP, T_POP, T_DO2) ;

# Because AIBECS will solve for the steady state solution directly without time-stepping the goverining equations to equilibrium, we don't have any opportunity to specify any intial conditions.
# Initial conditions are how the total amount of conserved elements get specified in most global biogeochemical modelels.
# Thus to specify the total inventory of P in AIBECS we add a very weak resporing term to the DIP equation.
# The time-scale for this restoring term is chosen to be very long compared to the timescale with which the ocean circulation homogenizes a tracer.
# Because of this long timescale we call it the geological restoring term, but geochemists who work on geological processes don't like that name!
# In any event the long timescale allows us to prescribe the total inventory of P in a way that yields the same solution we would have gotten had we time-stepped the model to steady-state with the total inventory prescribed by the initial condition.

# ### Sources minus sinks

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

function respiration(DOP, p)
    κDOP = p.κDOP
    return κDOP * DOP
end

# ##### Dissolution of POP into DOP

function dissolution(POP, p)
    κPOP = p.κPOP
    return κPOP * POP
end

# ##### Air-sea gas exchange

h = ustrip(grd.δdepth[1])             # thickness of the top layer
z = grd.depth_3D[iwet]        # depth of the gridbox centers
using WorldOceanAtlasTools
WOA = WorldOceanAtlasTools
μDO2sat , σ²DO2sat = WOA.fit_to_grid(grd,2018,"O2sat","annual","1°","an") ;
DO2sat = vec(μDO2sat)[iwet]
function airsea(DO2, p)
    κDO2 = p.κDO2
    return κDO2 * (z .< 20u"m") .* (DO2sat .- DO2) / h
end

# Add them up into sms functions (Sources Minus Sinks)

function sms_DIP(DIP, DOP, POP, DO2, p)
    return -uptake(DIP, p) + respiration(DOP, p) + geores(DIP, p)
end
function sms_DOP(DIP, DOP, POP, DO2, p)
    σ = p.σ
    return σ * uptake(DIP, p) - respiration(DOP, p) + dissolution(POP, p)
end
function sms_POP(DIP, DOP, POP, DO2, p)
    σ = p.σ
    return (1 - σ) * uptake(DIP, p) - dissolution(POP, p)
end
function sms_DO2(DIP, DOP, POP, DO2, p)
    rO2P = p.rO2P
    return airsea(DO2,p) + rO2P * (uptake(DIP,p) - respiration(DOP,p))
end
sms_all = (sms_DIP, sms_DOP, sms_POP, sms_DO2,) # bundles all the source-sink functions in a tuple


# ### Parameters
#
# Build the parameters type and p₀

t = empty_parameter_table()    # initialize table of parameters
add_parameter!(t, :xgeo, 2.17u"mmol/m^3")
add_parameter!(t, :τg, 1.0u"Myr")
add_parameter!(t, :ku, 10.0u"μmol/m^3")
add_parameter!(t, :z₀, 80.0u"m")
add_parameter!(t, :w₀, 1.0u"m/d")
add_parameter!(t, :w′, 1/4.4625u"d")
add_parameter!(t, :κDOP, 1/0.25u"yr")
add_parameter!(t, :κPOP, 1/5.25u"d")
add_parameter!(t, :σ, 0.3u"1")
add_parameter!(t, :τu, 30.0u"d")
add_parameter!(t, :κDO2, 50u"m" / 30u"d")
add_parameter!(t, :rO2P, 175.0u"mol/mol")
initialize_Parameters_type(t, "DIP_DOP_POP_O₂_Parameters")   # Generate the parameter type

# ### Generate state function and Jacobian

nt = length(T_all)    # number of tracers
n = nt * nb           # total dimension of the state vector
p = DIP_DOP_POP_O₂_Parameters() # parameters
x = p.xgeo * ones(n) # initial iterate
F, ∇ₓF = state_function_and_Jacobian(T_all, sms_all, nb)

# and solve

prob = SteadyStateProblem(F, ∇ₓF, x, p)
nothing # s = solve(prob, CTKAlg()) # Not working yet

# unpack state

DIP, DOP, POP, DO2 = state_to_tracers(x, nb, nt) # remove when line below works
nothing # DIP, DOP, POP, DO2 = state_to_tracers(s.u, nb, nt)

# We will plot the concentration of DIP at a given depth horizon

depth = grd.depth
iz = findfirst(depth .> 200u"m")
iz, depth[iz]

# Unfinished example! TODO

DIP_3D = rearrange_into_3Darray(DIP, grd)
DIP_2D = DIP_3D[:,:,iz] * ustrip(1.0u"mol/m^3" |> u"mmol/m^3")
lat, lon = ustrip.(grd.lat), ustrip.(grd.lon)

# and plot


