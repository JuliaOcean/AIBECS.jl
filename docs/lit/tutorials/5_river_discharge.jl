
#---------------------------------------------------------
# # [River discharge](@id river-discharge)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/5_river_discharge.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/5_river_discharge.ipynb)

# In this tutorial we will simulate a fictitious radioactive tracer that is injected into the ocean by the 200 largest rivers (by estimated discharge).
# The 200 major rivers dataset from [*Dai and Trenberth* (2002)](https://rda.ucar.edu/datasets/ds551.0/) is available from within the AIBECS.
# Once "born", our ficitious tracer decays with a parameter timescale $\tau$ as it flows through ocean basins.
#
# The 3D tracer equation is:
#
# $$\left[\frac{\partial}{\partial t} + \nabla \cdot (\boldsymbol{u} + \mathbf{K}\nabla)\right] x = s_\mathsf{rivers} - x / \tau$$
#
# where $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \nabla \right]$ represents the ocean circulation transport.
# ([Tracer transport operators are described in the documentation](@ref tracer-transport-operators).)
# The riverine source of the tracer is $s_\mathsf{rivers}$, and $x / \tau$ is the decay rate.

# In AIBECS, we must recast this equation in the generic form
#
# $$\left[\frac{\partial}{\partial t} + \mathbf{T}(\boldsymbol{p})\right] \boldsymbol{x} = \boldsymbol{G}(\boldsymbol{x},\boldsymbol{p}).$$

# We start by telling Julia that we want to use the AIBECS and the OCIM2 transport matrix for the ocean circulation.

using AIBECS
grd, T_OCIM2 = OCIM2.load()

# The transport is

T_radiorivers(p) = T_OCIM2

# For the radioactive decay, we simply use

function decay(x, p)
    @unpack τ = p
    return x / τ
end

# To build the river sources, we will load the geographic locations and discharge (in m³ s⁻¹) from the [*Dai and Trenberth* (2017) dataset](https://rda.ucar.edu/datasets/ds551.0/index.html#!description).

RIVERS = Rivers.load()

# This is an array of rivers, for which the type `River{T}` contains the river's name, lat–lon coordinates, and discharge in m³ s⁻¹.
# For example, the first river is the Amazon

r = RIVERS[1]

# We can check the locations with

using Plots
scatter([r.lon for r in RIVERS], [r.lat for r in RIVERS],
        zcolor=ustrip.([r.VFR for r in RIVERS] / u"m^3/s"),
        colorbartitle="log₁₀(discharge / (1 m³ s⁻¹))")

# We can regrid these into the OCIM2 grid and return the corresponding vector with

rivers = regrid(RIVERS, grd)

# (Note this regridding uses [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) to assign a wet box as the mouth of each river, which sometimes is not exactly the real location of the river mouth.)

# We control the global magnitude of the river discharge, $\sigma$ (in mol s⁻¹), by making it a parameter of our model.
# For that, we separate the river source
#
# $$s_\mathsf{rivers} = \sigma s_0$$
#
# into global magnitude ($\sigma$) and spatial pattern ($s_0$).
#
# Since $\int s_0 \mathrm{d}V = 1$, `s_0` can be computed by normalizing `rivers`.
# In Julia/AIBECS, this can be done by dividing `rivers` by the dot product `v ⋅ rivers` (or `v'rivers` in matrix form).
# (`v ⋅ x` is the discrete equivalent of the volume integral $\int x \mathrm{d}V$.)

v = vector_of_volumes(grd)
s_0 = rivers / (v'rivers)
function s_rivers(p)
    @unpack σ = p
    return σ * ustrip.(s_0) # we must remove the units in AIBECS here :(
end

# We then write the generic $\boldsymbol{G}$ function, which is

G_radiorivers(x,p) = s_rivers(p) - decay(x,p)

# ##### Parameters

# We specify some initial values for the parameters and also include units.

import AIBECS: @units, units
import AIBECS: @initial_value, initial_value
@initial_value @units struct RadioRiversParameters{U} <: AbstractParameters{U}
    τ::U | 5.0 | u"yr"
    σ::U | 1.0 | u"Gmol/yr"
end

# Finally, thanks to the initial values we provided, we can instantiate the parameter vector succinctly as

p = RadioRiversParameters()

# We generate the state function `F` and its Jacobian `∇ₓF`,

F, ∇ₓF = state_function_and_Jacobian(T_radiorivers, G_radiorivers)

# generate the steady-state problem `prob`,

nb = sum(iswet(grd))
x = ones(nb) # initial guess
prob = SteadyStateProblem(F, ∇ₓF, x, p)

# and solve it

s = solve(prob, CTKAlg()).u * u"mol/m^3"

# Let's now run some visualizations using the plot recipes.
# Taking a horizontal slice of the 3D field at 200m gives

cmap = :viridis
plothorizontalslice(s, grd, zunit=u"μmol/m^3", depth=200, color=cmap, clim=(0,0.5))

# and at 500m:

plothorizontalslice(s, grd, zunit=u"μmol/m^3", depth=500, color=cmap, clim=(0,0.05))

# Or we can change the timescale and watch the tracer fill the oceans:

p = RadioRiversParameters(τ = 50.0u"yr")
prob = SteadyStateProblem(F, ∇ₓF, x, p)
s_τ50 = solve(prob, CTKAlg()).u * u"mol/m^3"
plothorizontalslice(s_τ50, grd, zunit=u"μmol/m^3", depth=500, color=cmap, clim=(0,1))
