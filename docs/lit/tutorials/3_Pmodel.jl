
#---------------------------------------------------------
# # [A coupled PO₄–POP model](@id P-model)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/3_Pmodel.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/3_Pmodel.ipynb)

# In this tutorial we will explicitly simulate 2 tracers whose distributions control and feed back on each other.

# We consider a simple model for the cycling of phosphorus with 2 state variables consisting of phosphate (PO₄) AKA dissolved inorganic phosphorus (DIP) and particulate organic phosphorus (POP).
# The dissolved phases are transported by advection and diffusion whereas the particulate phase sinks rapidly down the water column without any appreciable transport by the circulation.
#
# The governing equations that couple the 3D concentration fields of DIP and POP, denoted $x_\mathsf{DIP}$ and $x_\mathsf{POP}$, respectively, are:
#
# $$\left[\frac{\partial}{\partial t} + \nabla \cdot (\boldsymbol{u} + \mathbf{K}\nabla )\right] x_\mathsf{DIP} = -U(x_\mathsf{DIP}) + R(x_\mathsf{POP}),$$
#
# and
#
# $$\left[\frac{\partial}{\partial t} + \nabla \cdot \boldsymbol{w}\right] x_\mathsf{POP} = U(x_\mathsf{DIP}) - R(x_\mathsf{POP}).$$

# The $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \nabla \right]$ and $\nabla \cdot \boldsymbol{w}$ operators represent the ocean circulation and the sinking of particles, respectively.
# ([Tracer transport operators are described in the documentation](@ref tracer-transport-operators).)

# The function $U$ represents the biological uptake of DIP by phytoplankton, which we model here as
#
# $$U(x_\mathsf{DIP}) = \frac{x_\mathsf{DIP}}{\tau_\mathsf{DIP}} \, \frac{x_\mathsf{DIP}}{x_\mathsf{DIP} + k} \, (z < z_0),$$
#
# with the timescale, $\tau$, the half-saturation rate $k$, and the depth $z_0$ as parameters.

# The function $R$ defines the remineralization rate of POP, which converts POP back into DIP.
# For the remineralization, we simply use a linear rate constant, i.e.,
#
# $$R(x_\mathsf{POP}) = \frac{x_\mathsf{POP}}{\tau_\mathsf{POP}}.$$



# We start by telling Julia we want to use the AIBECS and the OCIM0.1 circulation for DIP.

using AIBECS
grd, T_OCIM = OCIM0.load()
T_DIP(p) = T_OCIM

# For the sinking of particles, we use the `transportoperator` function

T_POP(p) = transportoperator(grd, z -> w(z,p))

# for which we need to define the sinking speed `w(z,p)` as a function of depth `z` and of the parameters `p`.
# Following the assumption that $w(z) = w_0 + w' z$ increases linearly with depth, we write it as

function w(z,p)
    @unpack w₀, w′ = p
    return @. w₀ + w′ * z
end

# ##### Uptake (DIP → POP)

# For the uptake, $U$, we write

z = depthvec(grd)
function U(x,p)
    @unpack τ_DIP, k, z₀ = p
    return @. x/τ_DIP * x/(x+k) * (z≤z₀) * (x≥0)
end

# where we have "unpacked" the parameters to make the code clearer and as close to the mathematical equation as possible.
# (Note we have also added a constraint that `x` must be positive for uptake to happen.)

# ##### Remineralization (POP → DIP)

# For the remineralization, $R$, we write

function R(x,p)
    @unpack τ_POP = p
    return x / τ_POP
end



# ##### Net sources and sinks

# We lump the sources and sinks into `G` functions for DIP and POP.

function G_DIP(DIP, POP, p)
    @unpack DIP_geo, τ_geo = p
    return @. -$U(DIP,p) + $R(POP,p) + (DIP_geo - DIP) / τ_geo
end
function G_POP(DIP, POP, p)
    @unpack τ_geo = p
    return @. $U(DIP,p) - $R(POP,p) - POP / τ_geo
end

# where we have imposed a slow restoring of DIP to the global mean `DIP_geo` to prescribe the global mean concentration.
# (The `$` signs in front of `U` and `R` protect them from the broadcast macro `@.`)

# We now define and build the parameters.

# In this tutorial we will specify some initial values for the parameters
# and also include units.

import AIBECS: @units, units
import AIBECS: @initial_value, initial_value
using Unitful: m, d, s, yr, Myr, mol, mmol, μmol, μM
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

# Finally, thanks to the initial values we provided, we can instantiate the parameter vector succintly as

p = PmodelParameters()

# We generate the state function `F` and its Jacobian `∇ₓF`,

nb = sum(iswet(grd))
F, ∇ₓF = state_function_and_Jacobian((T_DIP, T_POP), (G_DIP, G_POP), nb)

# generate the steady-state problem,

@unpack DIP_geo = p
x = DIP_geo * ones(2nb) # initial guess
prob = SteadyStateProblem(F, ∇ₓF, x, p)

# and solve it

sol = solve(prob, CTKAlg()).u

# We can look at different the DIP and POP fields using the Plots.jl recipes.

DIP, POP = state_to_tracers(sol, grd) # unpack tracers

# First, let's look at the mean profile

using Plots
plothorizontalmean(DIP * (mol/m^3) .|> μM, grd)

# We can plot the concentration of DIP at a given depth via, e.g.,

plothorizontalslice(DIP * (mol/m^3) .|> μM, grd, depth=1000m, color=:viridis)

# Or have a look at a map of the uptake at the surface

plotverticalintegral(U(DIP,p) * (mol/m^3/s) .|> mmol/yr/m^3, grd, color=:algae)

# Or look at what is exported below 500 m

plothorizontalslice(POP .* w(z,p) * (mol/m^3*m/s) .|> mmol/yr/m^2, grd, depth=500m, color=:inferno, rev=true)

# Now let's make our model a little fancier and use a fine topographic map to refine the remineralization profile.
# For this, we will use the ETOPO dataset, which can be downloaded by AIBECS via

f_topo = ETOPO.fractiontopo(grd)

# `f_topo` is the fraction of sediments located in each wet box of the `grd` grid.
# We can use it to redefine the transport operator for sinking particles to take into consideration the subgrid topography,
# such that the fine-resolution sediments intercept settling POP.

T_POP2(p) = transportoperator(grd, z -> w(z,p); frac_seafloor=f_topo)

# With this new vertical transport for POP, we can recreate our problem, solve it again

F2, ∇ₓF2 = state_function_and_Jacobian((T_DIP, T_POP2), (G_DIP, G_POP), nb)
prob2 = SteadyStateProblem(F2, ∇ₓF2, x, p)
sol2 = solve(prob2, CTKAlg()).u
DIP2, POP2 = state_to_tracers(sol2, grd) # unpack tracers

# and check the difference

plotzonalaverage((DIP2 - DIP) ./ DIP .|> u"percent", grd, color=:balance, clim=(-0.1, 0.1))

# This zonal average shows how much DIP is prevented from sinking out of the surface layers with the new subgrid parameterization.

# Let's look at the vertical average.

plotverticalaverage((DIP2 - DIP) ./ DIP .|> u"percent", grd, color=:balance, clim=(-0.1,0.1))

# This shows minor changes, on the order of 0.1%, on the global scale,
# which means that the subgrid-topography parameterization mostly affects the vertical distribution of our tracers.
