
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

T_POP(p) = transportoperator(grd, w = w(p))

# for which we need to define the sinking speed `w(p)` as a function of the parameters `p`.
# Following the assumption that $w = w_0 + w' z$ increases linearly with depth, we write it as

function w(p)
    @unpack w₀, w′ = p
    return @. w₀ + w′ * z
end

# For this to work, we must create a vector of depths, `z`, which is simply done via

iwet = findall(vec(iswet(grd)))
z = ustrip.(grd.depth_3D[iwet])

# (We strip the units of `z` via `ustrip`.)


# ##### Uptake (DIP → POP)

# For the uptake, $U$, we write

function U(x,p)
    @unpack τDIP, k, z₀ = p
    return @. x/τDIP * x/(x+k) * (z≤z₀) * (x≥0)
end

# where we have "unpacked" the parameters to make the code clearer and as close to the mathematical equation as possible.
# (Note we have also added a constraint that `x` must be positive for uptake to happen.)

# ##### Remineralization (POP → DIP)

# For the remineralization, $R$, we write

function R(x,p)
    @unpack τPOP = p
    return x / τPOP
end



# ##### Net sources and sinks

# We lump the sources and sinks into `G` functions for DIP and POP.

function G_DIP(DIP, POP, p)
    @unpack D̅I̅P̅, τgeo = p
    return @. -$U(DIP,p) + $R(POP,p) + (D̅I̅P̅ - DIP) / τgeo
end
function G_POP(DIP, POP, p)
    return U(DIP,p) - R(POP,p)
end

# where we have imposed a slow restoring of DIP to the global mean `D̅I̅P̅` to prescribe the global mean concentration.
# (The `$` signs in front of `U` and `R` protect them from the broadcast macro `@.`)

# We now define and build the parameters.

# In this tutorial we will specify some initial values for the parameters
# and also include units.

import AIBECS: @units, units
import AIBECS: @initial_value, initial_value
@units @initial_value struct PmodelParameters{U} <: AbstractParameters{U}
    w₀::U   |  0.64 | u"m/d"
    w′::U   |  0.13 | u"m/d/m"
    τDIP::U | 230.0 | u"d"
    k::U    |  6.62 | u"μmol/m^3"
    z₀::U   |  80.0 | u"m"
    τPOP::U |   5.0 | u"d"
    τgeo::U |   1.0 | u"Myr"
    D̅I̅P̅::U  |  2.12 | u"mmol/m^3"
end

# Finally, thanks to the initial values we provided, we can instantiate the parameter vector succintly as

p = PmodelParameters()

# We generate the state function `F` and its Jacobian `∇ₓF`,

nb = sum(iswet(grd))
F, ∇ₓF = state_function_and_Jacobian((T_DIP, T_POP), (G_DIP, G_POP), nb)

# generate the steady-state problem,

@unpack D̅I̅P̅ = p
x = D̅I̅P̅ * ones(2nb) # initial guess
prob = SteadyStateProblem(F, ∇ₓF, x, p)

# and solve it

s = solve(prob, CTKAlg()).u

# We can look at different the DIP and POP fields using the Plots.jl recipes.

DIP, POP = state_to_tracers(s, grd) # unpack tracers

# We can plot the concentration of DIP at a given depth via, e.g.,

using Plots
horizontalslice(DIP * u"mol/m^3" .|> u"μM", grd, 1000; color=:viridis)

# Or have a look at a map of the uptake at the surface

verticalintegral(U(DIP,p) * u"mol/m^3/s" .|> u"mmol/yr/m^3", grd; color=:algae)

# Or look at what is exported below 500 m

horizontalslice(POP .* w(p) * u"mol/m^3*m/s" .|> u"mmol/yr/m^2", grd, 500; color=:inferno_r)

