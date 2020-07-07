
#---------------------------------------------------------
# # [A dust model](@id river-discharge)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/5_river_discharge.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/5_river_discharge.ipynb)

# In this tutorial we will simulate a fictitious radioactive tracer that is injected into the ocean by the 200 largest rivers (by estaimated discharge).
# The 200 major rivers are available from within the AIBECS.
# Once "born", our tracer will decay with a parameter timescale $\tau$ as it flows through ocean basins.
#
# The 3D tracer equation is:
#
# $$\left[\frac{\partial}{\partial t} + \nabla \cdot (\boldsymbol{u} + \mathbf{K}\nabla)\right] x = s_\mathsf{rivers} - x / \tau$$
#
# where $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \nabla \right]$ represents the ocean circulation transport.
# ([Tracer transport operators are described in the documentation](@ref tracer-transport-operators).)
# The riverine source of the tracer is $s_\mathsf{rivers}$, and it decays with rate $x / \tau$.

# In AIBECS, we must recast this equation in the generic form
#
# $$\left[\frac{\partial}{\partial t} + \mathbf{T}(\boldsymbol{p})\right] \boldsymbol{x} = \boldsymbol{G}(\boldsymbol{x},\boldsymbol{p}).$$

# We start by telling Julia we want to use the AIBECS and the OCIM0.1 transport matrix for the ocean circulation.

using AIBECS
grd, T_OCIM = OCIM0.load()

# The transport is

T_radiorivers(p) = T_OCIM

# For the radioactive decay, we simply use

function decay(x, p)
    @unpack τ = p
    return x / τ
end

# To build the river sources, we will load the geographic locations and discharge (in m<sup>3</sup> s<sup>-1</sup>) from the [Dai and Trenberth (2017) dataset](https://rda.ucar.edu/datasets/ds551.0/index.html#!description).

RIVERS = Rivers.load()

# This is an array of rivers, for which the type`River{T}` contains the river's name, lat–lon coordinates, an discharge in m<sup>3</sup> s<sup>-1</sup>.

# We can regrid these into the OCIM0.1 grid and return the corresponding vector with

rivers = regrid(RIVERS, grd)

# (Note this regridding uses [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) to assigne a wet box as the mouth of each river, which sometimes is not exactly the real loaction of the mouth.)
# We then normalize the river discharge globally to control the global magnitude of the riverine sources with a single parameter, $\sigma$ (in mol s<sup>-1</sup>.
# We want $\int s_\mathsf{rivers} dV = \sigma_\mathsf{rivers}$ (in mol s<sup>-1</sup>).
# We want to write it as $s_\mathsf{rivers} = \sigma_\mathsf{rivers} s_\mathsf{rivers}'$, so $\int s_\mathsf{rivers}' dV = 1$. In Julia/AIBECS, this is equivalent to `v's_rivers′ == 1` where `v` is the vector of volumes. Thus, we do

v = vector_of_volumes(grd)
s_rivers′ = rivers / (v'rivers)

# assuming
# and scale
function s_rivers(p)
    @unpack σ = p
    return σ * ustrip.(s_rivers′) # we must remove the units in AIBECS here :(
end

# We then write the generic $\boldsymbol{G}$ function, which is

G_radiorivers(x,p) = s_rivers(p) - decay(x,p)

# ##### Parameters

# We specify some initial values for the parameters and also include units.

import AIBECS: @units, units
import AIBECS: @initial_value, initial_value
@initial_value @units struct RadioRiversParameters{U} <: AbstractParameters{U}
    τ::U        | 5.0 | u"yr"
    σ::U        | 1.0 | u"Gmol/yr"
end

# Finally, thanks to the initial values we provided, we can instantiate the parameter vector succintly as

p = RadioRiversParameters()

# We generate the state function `F` and its Jacobian `∇ₓF`,

F, ∇ₓF = state_function_and_Jacobian(T_radiorivers, G_radiorivers)

# generate the steady-state problem,

nb = sum(iswet(grd))
x = ones(nb) # initial guess
prob = SteadyStateProblem(F, ∇ₓF, x, p)

# and solve it

s = solve(prob, CTKAlg()).u * u"mol/m^3"

# Let's now run some vizualizations using the Plots.jl recipes.

using Plots

# Let's have a look at the distribution layer by layer.
# At 100m,

plothorizontalslice(s, grd, zunit=u"μmol/m^3", depth=100, color=:linear_bmw_5_95_c89_n256)

# and 500m.

plothorizontalslice(s, grd, zunit=u"μmol/m^3", depth=500, color=:linear_bmw_5_95_c89_n256)

# Or we can change the timescale and watch the tracer fill the oceans.

p = RadioRiversParameters(τ = 50.0u"yr")
prob = SteadyStateProblem(F, ∇ₓF, x, p)
s_τ50 = solve(prob, CTKAlg()).u * u"mol/m^3"
plothorizontalslice(s_τ50, grd, zunit=u"μmol/m^3", depth=500, color=:linear_bmw_5_95_c89_n256)
