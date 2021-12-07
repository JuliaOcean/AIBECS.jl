
#---------------------------------------------------------
# # [Radiocarbon](@id radiocarbon)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/2_radiocarbon.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/2_radiocarbon.ipynb)

#md # !!! note
#md #     All the AIBECS tutorials and how-to guides are available as Jupyter notebooks.
#md #     You can view them with [nbviewer](https://nbviewer.jupyter.org/)or execute them online with [binder](https://mybinder.org/) by clicking on the badges above!
#md #     (Note that binder can be slow to launch and its memory caps can be a problem when running.)

# In this tutorial, we will simulate the radiocarbon age using the AIBECS by
# 1. defining the transport `T(p)` and the sources and sinks `G(x,p)`,
# 1. defining the parameters `p`,
# 1. generating the state function `F(x,p)` and solving the associated steady-state problem,
# 1. and finally making a plot of our simulated radiocarbon age.


#md # !!! note
#md #     Although this tutorial is self-contained, it is slightly more complicated than [the first tutorial for simulating the ideal age](@ref ideal-age).
#md #     (So do not hesitate to start with the ideal-age tutorial if you wish.)

#nb # > *Note*
#nb # > Although this tutorial is self-contained, it is slightly more complicated than [the first tutorial for simulating the ideal age](@ref ideal-age).
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

# We start by selecting the circulation for Radiocarbon.
# .)
#nb # (And this time, we are using the OCCA matrix by *Forget* [1](https://doi.org/10.1175/2009JPO4043.1).)
#md # (And this time, we are using the OCCA matrix by *Forget* [^1].)

#md # [^1]:
#md #     Forget, G., 2010: Mapping Ocean Observations in a Dynamical Framework: A 2004–06 Ocean Atlas. J. Phys. Oceanogr., 40, 1201–1221, doi:[10.1175/2009JPO4043.1)](https://doi.org/10.1175/2009JPO4043.1)


using AIBECS
grd, T_OCCA = OCCA.load()

# The local sources and sinks are simply given by

function G(R,p)
    @unpack λ, h, Ratm, τ = p
    return @. λ / h * (Ratm - R) * (z ≤ h) - R / τ
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
    Ratm::U | u"M"
end

# For the air–sea gas exchange, we use a constant piston velocity $\lambda$ of 50m / 10years.
# And for the radioactive decay we use a timescale $\tau$ of 5730/log(2) years.
#

p = RadiocarbonParameters(λ = 50u"m"/10u"yr",
                          h = grd.δdepth[1],
                          τ = 5730u"yr"/log(2),
                          Ratm = 42.0u"nM")

#md # !!! note
#md #     The parameters are converted to SI units when unpacked.
#md #     When you specify units for your parameters, you must either supply their values in that unit.
#nb # > *Note*
#nb # > The parameters are converted to SI units when unpacked.
#nb # >  When you specify units for your parameters, you must supply their values in that unit.

# We build the state function `F` and the corresponding steady-state problem (and solve it) via

F = AIBECSFunction(T_OCCA, G)
x = zeros(length(z)) # an initial guess
prob = SteadyStateProblem(F, x, p)
R = solve(prob, CTKAlg()).u

# This should take a few seconds on a laptop.
# Once the radiocarbon concentration is computed, we can convert it into the corresponding age in years, via

@unpack τ, Ratm = p
C14age = @. log(Ratm / R) * τ * u"s" |> u"yr"

# and plot it at 700 m using the `horizontalslice` Plots recipe

using Plots
plothorizontalslice(C14age, grd, depth=700u"m", color=:viridis)

# look at a zonal average using the `zonalaverage` plot recipe

plotzonalaverage(C14age, grd; color=:viridis)

# or look at a meridional slice through the Atlantic at 30°W using the `meridionalslice` plot recipe

plotmeridionalslice(C14age, grd, lon=-30, color=:viridis)
