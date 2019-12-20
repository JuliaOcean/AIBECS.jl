
#---------------------------------------------------------
# # [Radiocarbon](@id radiocarbon)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/2_radiocarbon.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/2_radiocarbon.ipynb)

# In this tutorial, we will simulate the radiocarbon age using the AIBECS by
# 1. defining the transport `T(p)` and the sources and sinks `G(x,p)`,
# 1. defining the parameters `p`,
# 1. generating the state function `F(x,p)` and solving the associated steady-state problem,
# 1. and finally making a plot of our simulated radiocarbon age.


#md # !!! note
#md #     Although this tutorial is self-contained, it involves non-linearitiess and is slightly more complicated than [the first tutorial for simulating the ideal age](@ref ideal-age).
#md #     (So do not hesitate to start with the ideal-age tutorial if you wish.)

#nb # > *Note*
#nb # > Although this tutorial is self-contained, it involves non-linearitiess and is slightly more complicated than [the first tutorial for simulating the ideal age](@ref ideal-age).
#nb # > (So do not hesitate to start with the ideal-age tutorial if you wish.)

# The tracer equation for radiocarbon is
#
# $$\big(\partial_t + \mathbf{T} \big) \boldsymbol{R} = \frac{\lambda}{h} (\overline{\boldsymbol{R}}_\mathsf{atm} - \boldsymbol{R}) (\boldsymbol{z} ≤ h) - \boldsymbol{R} / \tau.$$
#
# where the first term on the right of the equal sign represents the air–sea gas exchange with a piston velocity $λ$ over a depth $h$ and the second term represents the radioactive decay of radiocarbon with timescale $\tau$.



#md # !!! note
#md #     We need not specify the value of the atmospheric radiocarbon concentration because it is not important for determining the age of a water parcel — only the relative concentration $\boldsymbol{R}/\overline{\boldsymbol{R}}_\mathsf{atm}$ matters.

#nb # > Note:
#nb # > We need not specify the value of the atmospheric radiocarbon concentration because it is not important for determining the age of a water parcel — only the relative concentration $\boldsymbol{R}/\overline{\boldsymbol{R}}_\mathsf{atm}$ matters.

# We start by selecting the circulation for Radiocarbon

using AIBECS
grd, T_OCIM1 = OCIM1.load()
T(p) = T_OCIM1

# The local sources and sinks are simply given by

function G(R,p)
    @unpack λ, h, R̅atm, τ = p
    return @. λ / h * (R̅atm - R) * (z ≤ h) - R / τ
end


# We can define `z` via

z = depthvec(grd)

# In this tutorial we will specify some units for the parameters.
# Such features **must** be imported to be used

import AIBECS: @units, units

# We define the parameters using the dedicated API from the AIBECS, including keyword arguments and units this time

@units struct RadiocarbonParameters{U} <: AbstractParameters{U}
    λ::U    | u"m/yr"
    h::U    | u"m"
    τ::U    | u"yr"
    R̅atm::U | u"M"
end

# For the air–sea gas exchange, we use a constant piston velocity $\lambda$ of 50m / 10years.
# And for the radioactive decay we use a timescale $\tau$ of 5730/log(2) years.
#

p = RadiocarbonParameters(λ = 50u"m"/10u"yr", 
                          h = grd.δdepth[1], 
                          τ = 5730u"yr"/log(2), 
                          R̅atm = 42.0u"nM")

#md # !!! note
#md #     The parameters are converted to SI units when unpacked.
#md #     When you specify units for your parameters, you must either supply their values in that unit.
#nb # > *Note*
#nb # > The parameters are converted to SI units when unpacked.
#nb # >  When you specify units for your parameters, you must supply their values in that unit.

# We generate the state function and its Jacobian, generate the corresponding steady-state problem, and solve it, via

F, ∇ₓF = state_function_and_Jacobian(T, G)
x = zeros(length(z)) # an initial guess
prob = SteadyStateProblem(F, ∇ₓF, x, p)
R = solve(prob, CTKAlg()).u

# This should take a few seconds on a laptop.
# Once the radiocarbon concentration is computed, we can convert it into the corresponding age in years, via

@unpack τ, R̅atm = p
C14age = @. log(R̅atm / R) * τ * u"s" |> u"yr"

# and plot it at 700 m using the `horizontalslice` Plots recipe

using Plots
horizontalslice(C14age, grd, 700; color=:viridis)

# look at a zonal average using the `zonalaverage` plot recipe

zonalaverage(C14age, grd; color=:viridis)

# or look at a zonal slice through the Atlantic at 30°W using the `zonalslice` plot recipe

zonalslice(C14age, grd, -30; color=:viridis)
