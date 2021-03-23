
#---------------------------------------------------------
# # [The pseudo-transient continuation solver](@id PTCsolver)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/7_PTCsolver.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/7_PTCsolver.ipynb)

# In this tutorial we use the pseudo-transient continuation solver to find the steady-state solution of a stiff 2-tracer model that is hard to solve directly with Newton's method.

# TODO: Fix text
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

# TODO: Fix text
# The $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \nabla \right]$ and $\nabla \cdot \boldsymbol{w}$ operators represent the ocean circulation and the sinking of particles, respectively.
# ([Tracer transport operators are described in the documentation](@ref tracer-transport-operators).)

# The function $U$ represents the biological uptake of DIP by phytoplankton, which we model here as
#
# $$U(x_\mathsf{DIP}) = \frac{x_\mathsf{DIP}}{\tau_\mathsf{DIP}} \, \frac{x_\mathsf{DIP}}{x_\mathsf{DIP} + k} \, (z < z_0),$$
#
# with the timescale, $\tau$, the half-saturation rate $k$, and the depth $z_0$ as parameters.

# TODO: Fix text
# The function $R$ defines the remineralization rate of POP, which converts POP back into DIP.
# For the remineralization, we simply use a linear rate constant, i.e.,
#
# $$R(x_\mathsf{POP}) = \frac{x_\mathsf{POP}}{\tau_\mathsf{POP}}.$$



# We start by telling Julia we want to use the AIBECS and the OCCA circulation for DIP.

using AIBECS
grd, T_OCCA = OCCA.load()
T_DIP = T_OCCA

# For the sinking of particles, we use the `transportoperator` function

function w(p)
    @unpack w₀, w′ = p
    z -> w₀ + w′ * z
end
T_POP = T_OCCA + transportoperator(grd, w)

# for which we need to define the sinking speed `w(z,p)` as a function of depth `z` and of the parameters `p`.
# Following the assumption that $w(z) = w_0 + w' z$ increases linearly with depth, we write it as


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
F, ∇ₓF = F_and_∇ₓF((T_DIP, T_POP), (G_DIP, G_POP), nb)

# generate the steady-state problem,

@unpack DIP_geo = p
x = DIP_geo * ones(2nb) # initial guess
prob = SteadyStateProblem(F, ∇ₓF, x, p)

# and solve it

# TODO use PTC solver instead here
sol = solve(prob, CTKAlg()).u






# 


