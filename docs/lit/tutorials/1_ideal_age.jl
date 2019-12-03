
#---------------------------------------------------------
# # [Ideal age](@id ideal-age)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/1_ideal_age.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/1_ideal_age.ipynb)

# The tracer equation for the ideal age is
#
# $$\left(\partial_t + \mathbf{T}\right) \boldsymbol{a} = 1 - \frac{\boldsymbol{a}}{τ} \, (\boldsymbol{z} < z_0),$$
#
# where the sink term on the right clamps the age to $0$ at the surface (where $\boldsymbol{z} < z_0$).
# The smaller the timescale $\tau$, the quicker $\boldsymbol{a}$ is restored to $0$ at the surface.

# AIBECS can interpret tracer equations as long as you arrange them under the generic form:
#
# $$\big(\partial_t + \mathbf{T}(\boldsymbol{p}) \big) \boldsymbol{x} = \boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p}),$$
#
# where $\mathbf{T}(\boldsymbol{p})$ is the transport, $\boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p})$ is the net local sources and sinks, and $\boldsymbol{p}$ is the vector of model parameters.
# We will then use the AIBECS to simulate the ideal age by finding the steady-state of the system, i.e., the solution of
#
# $$\partial_t \boldsymbol{x} = \boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = \boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p}) - \mathbf{T}(\boldsymbol{p}) \, \boldsymbol{x} = 0.$$

# In this tutorial, we will simulate the ideal age using the AIBECS by
# 1. defining functions for `T(p)` and `G(x,p)`,
# 1. defining the parameters `p`,
# 1. then generating the state function `F(x,p)` and solving the associated steady-state problem,
# 1. and finally making a plot of our simulated ideal age.

# We start by telling Julia that we want to use the AIBECS package and the OCIM1 circulation
#nb # (the Ocean Circulation Inverse Model, see [1](https://doi.org/10.1002/2013GB004739) and [2](https://doi.org/10.1175/JPO-D-10-05011.1) for details).
#md # (the Ocean Circulation Inverse Model[^1][^2]).

#md # [^1]:
#md #     DeVries, T. (2014), The oceanic anthropogenic CO₂ sink: Storage, air–sea fluxes, and transports over the industrial era, Global Biogeochem. Cycles, 28, 631–647, doi:[10.1002/2013GB004739](https://doi.org/10.1002/2013GB004739).

#md # [^2]:
#md #     DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, doi:[10.1175/JPO-D-10-05011.1](https://doi.org/10.1175/JPO-D-10-05011.1)

using AIBECS
grd, TOCIM1 = OCIM1.load()

#md # !!! note
#md #     If it's the first time, Julia may ask you to download the OCIM1, in which case you should accept (i.e., type `y` and "return").
#md #     Once downloaded, AIBECS will remember where it downloaded the file and it will only load it from your laptop.
#nb # > **Note**
#nb # > If it's the first time, Julia may ask you to download the OCIM1, in which case you should accept (i.e., type `y` and "return").
#nb # > Once downloaded, AIBECS will remember where it downloaded the file and it will only load it from your laptop.

# `grd` is an `OceanGrid` object containing information about the 3D grid of the OCIM1 circulation and `TOCIM` is the transport matrix representing advection and diffusion.

# We define the function `T(p)` as

T(p) = TOCIM1

# which turns out to not effectively depend on `p`, but that's how we must define `T(p)` anyway.

# The local sources and sinks for the age take the form

function G(x,p)
    @unpack τ, z₀ = p
    return @. 1 - x / τ * (z < z₀)
end

# as per the tracer equation.
# The `@unpack` line unpacks the parameters `τ` and `z₀`.
# The `return` line returns the net sources and sinks.
# (The `@.` "macro" tells Julia that the operations apply to every element.)

# We need to define the vector `z` of depths.

z = ustrip.(vector_of_depths(grd))

# (We stripped `z` from its units (m) with `ustrip`.)

# Now we must construct a type for `p` the parameters.
# This type must contain our parameters `τ` and `z₀`.

struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end

# The type is now ready for us to generate an instance of the parameter `p`.
# Let's use `τ = 1.0` (s) and `z₀ = 20.0` (m) .

p = IdealAgeParameters(1.0, 20.0)

# We now use the AIBECS to generate the state function $\boldsymbol{F}$ (and its Jacobian) via

F, ∇ₓF = state_function_and_Jacobian(T, G)

# (`∇ₓF` is the **Jacobian** of the state function $\nabla_{\boldsymbol{x}}\boldsymbol{F}$, calculated automatically using dual numbers.)

# Now that `F(x,p)`, and `p` are defined, we are going to solve for the steady-state.
# But first, we must create a `SteadyStateProblem` object that contains `F`, `∇ₓF`, `p`, and an initial guess `x_init` for the age.
# (`SteadyStateProblem` is specialized from [DiffEqBase](https://github.com/JuliaDiffEq/DiffEqBase.jl) for AIBECS models.)

# Let's make a vector of 0's for our initial guess.

nb = sum(iswet(grd))  # number of wet boxes
x_init = zeros(nb)    # Start with age = 0 everywhere

# Now we can create our `SteadyStateProblem` instance

prob = SteadyStateProblem(F, ∇ₓF, x_init, p)

# And finally, we can `solve` this problem, using the AIBECS `CTKAlg()` algorithm,

age = solve(prob, CTKAlg())

# This should take a few seconds.

# To conclude this tutorial, let's have a look at the age using AIBECS' plotting recipes and [Plots.jl](https://github.com/JuliaPlots/Plots.jl).

using Plots

# We first convert the age in years
# (because the default SI unit we used, i.e., seconds, is a bit small relative to global ocean timescales).

age_in_yrs = age * u"s" .|> u"yr"

# And we take a horizontal slice at about 2000m.

horizontalslice(age_in_yrs, grd, 2000; color=:magma)

# That's it for this tutorial...
# Good job!

