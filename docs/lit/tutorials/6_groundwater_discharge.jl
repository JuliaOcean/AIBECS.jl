
#---------------------------------------------------------
# # [Groundwater discharge](@id groundwater-discharge)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/6_groundwater_discharge.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/6_groundwater_discharge.ipynb)

#md # !!! note
#md #     This tutorial is very similar to the [river discharge notebook](@id river-discharge)

# In this tutorial we will simulate a fictitious radioactive tracer that is injected into the ocean by groundwater discharge.
# The global groundwater discharge dataset for 40,000 coastal watersheds from [*Luijendijk et al.* (2020)](https://www.nature.com/articles/s41467-020-15064-8) is available from within the AIBECS.
# Once "born", our ficitious tracer decays with a parameter timescale $\tau$ as it flows through ocean basins.
#
# The 3D tracer equation is:
#
# $$\left[\frac{\partial}{\partial t} + \nabla \cdot (\boldsymbol{u} + \mathbf{K}\nabla)\right] x = s_\mathsf{gw} - x / \tau$$
#
# where $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \nabla \right]$ represents the ocean circulation transport.
# ([Tracer transport operators are described in the documentation](@ref tracer-transport-operators).)
# The source of the tracer is $s_\mathsf{gw}$, and $x / \tau$ is the decay rate.

# In AIBECS, we must recast this equation in the generic form
#
# $$\left[\frac{\partial}{\partial t} + \mathbf{T}(\boldsymbol{p})\right] \boldsymbol{x} = \boldsymbol{G}(\boldsymbol{x},\boldsymbol{p}).$$

# We start by telling Julia that we want to use the AIBECS and the OCIM2 transport matrix for the ocean circulation.

using AIBECS
grd, T_OCIM2 = OCIM2.load()

# The transport is

T_gw(p) = T_OCIM2

# For the radioactive decay, we simply use

function decay(x, p)
    @unpack τ = p
    return x / τ
end

# To build the groundwater sources, we will load the geographic locations and discharge of groundwater (in m³ yr⁻¹) from the *Luijendijk et al.* (2020)) dataset.

gws = GroundWaters.load()

# This is an array of groundwater sources, for which the type `GroundWaterSource{T}` contains the lat–lon coordinates and discharge in m³ yr⁻¹.
# For example, the first element is

gws[1]

# We can check the locations with

using Plots
scatter([x.lon for x in gws], [x.lat for x in gws],
        zcolor=log10.(ustrip.([x.VFR for x in gws] / u"m^3/s")),
        colorbartitle="log₁₀(discharge / (1 m³ s⁻¹))",
        clim=(0,10), fmt=:png)

# We can regrid these into the OCIM2 grid and return the corresponding vector with

gws_grd = ustrip.(u"m^3/s", regrid(gws, grd))

# (We convert years to seconds for AIBECS)

# We can control the global magnitude of the discharge by chosing a concentration of our tracer in groundwater, $C_\mathsf{gw}$, in mol m⁻³.
# For that, we need the volumes of the grid boxes

v = volumevec(grd)

# because the source is given by the mass (discharge × concentration) divided by the volume of each box:

function s_gw(p)
    @unpack C_gw = p
    return @. gws_grd * C_gw / v
end

# We then write the generic $\boldsymbol{G}$ function, which is

G_gw(x,p) = s_gw(p) - decay(x,p)

# ##### Parameters

# We specify some initial values for the parameters and also include units.

import AIBECS: @units, units
import AIBECS: @initial_value, initial_value
@initial_value @units struct GroundWatersParameters{U} <: AbstractParameters{U}
    τ::U    | 20.0 | u"yr"
    C_gw::U |  1.0 | u"mol/m^3"
end

# Finally, thanks to the initial values we provided, we can instantiate the parameter vector succinctly as

p = GroundWatersParameters()

# We generate the state function `F` and its Jacobian `∇ₓF`,

F, ∇ₓF = state_function_and_Jacobian(T_gw, G_gw)

# generate the steady-state problem `prob`,

nb = sum(iswet(grd))
x = ones(nb) # initial guess
prob = SteadyStateProblem(F, ∇ₓF, x, p)

# and solve it

s = solve(prob, CTKAlg()).u * u"mol/m^3"

# Let's now run some visualizations using the plot recipes.
# Taking a horizontal slice of the 3D field at 200m gives

cmap = :viridis
plothorizontalslice(s, grd, zunit=u"μmol/m^3", depth=200, color=cmap, clim=(0,100))

# and at 500m:

plothorizontalslice(s, grd, zunit=u"μmol/m^3", depth=500, color=cmap, clim=(0,25))

# Or we can increase the decay timescale (×10) and decrease the groundwater concentration (÷10) to get a different (more well-mixed) tracer distribution:

p = GroundWatersParameters(τ = 200.0u"yr", C_gw = 0.1u"mol/m^3")
prob = SteadyStateProblem(F, ∇ₓF, x, p)
s_τ50 = solve(prob, CTKAlg()).u * u"mol/m^3"
plothorizontalslice(s_τ50, grd, zunit=u"μmol/m^3", depth=500, color=cmap, clim=(0,25))
