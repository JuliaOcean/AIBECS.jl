
#---------------------------------------------------------
# # [Ideal age](@id ideal-age)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/1_ideal_age.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/1_ideal_age.ipynb)

#md # !!! note
#md #     All the AIBECS tutorials and how-to guides are available as Jupyter notebooks.
#md #     You can execute them online with [binder](https://mybinder.org/) or just view them with [nbviewer](https://nbviewer.jupyter.org/) by clicking on the badges above!

# The tracer equation for the ideal age is
#
# $$\left(\partial_t + \mathbf{T}\right) \boldsymbol{a} = 1 - \frac{\boldsymbol{a}}{τ} \, (\boldsymbol{z} \le z_0),$$
#
# where the sink term on the right clamps the age to $0$ at the surface (where $\boldsymbol{z} \le z_0$).
# The smaller the timescale $\tau$, the quicker $\boldsymbol{a}$ is restored to $0$ at the surface.

# AIBECS can interpret tracer equations as long as you arrange them under the generic form:
#
# $$\big(\partial_t + \mathbf{T}(\boldsymbol{p}) \big) \boldsymbol{x} = \boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p}),$$
#
# where $\mathbf{T}(\boldsymbol{p})$ is the transport, $\boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p})$ is the net local sources and sinks, and $\boldsymbol{p}$ is the vector of model parameters.
# We will then use the AIBECS to simulate the ideal age by finding the steady-state of the system, i.e., the solution of
#
# $$\partial_t \boldsymbol{x} = \boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = \boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p}) - \mathbf{T}(\boldsymbol{p}) \, \boldsymbol{x} = 0.$$

# In this tutorial, we will simulate the ideal age by
# 1. defining functions for `T(p)` and `G(x,p)`,
# 1. defining the parameters `p`,
# 1. generating the state function `F(x,p)` and solving the associated steady-state problem,
# 1. and finally making a plot of our simulated ideal age.

# We start by telling Julia that we want to use the AIBECS package and the OCIM2 circulation
#nb # (the Ocean Circulation Inverse Model, see [1](https://doi.org/10.1029/2018JC014716) for details).
#md # (the Ocean Circulation Inverse Model[^1]).

#md # [^1]:
#md #     DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle‐³He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. doi:[10.1029/2018JC014716](https://doi.org/10.1029/2018JC014716)

using AIBECS
grd, TOCIM2 = OCIM2.load()

#md # !!! note
#md #     If it's your first time, Julia will ask you to download the OCIM2, in which case you should accept (i.e., type `y` and "return").
#md #     Once downloaded, AIBECS will remember where it downloaded the file and it will only load it from your laptop.
#nb # > **Note**
#nb # > If it's your first time, Julia will ask you to download the OCIM2, in which case you should accept (i.e., type `y` and "return").
#nb # > Once downloaded, AIBECS will remember where it downloaded the file and it will only load it from your laptop.

# `grd` is an `OceanGrid` object containing information about the 3D grid of the OCIM2 circulation and `TOCIM2` is the transport matrix representing advection and diffusion.

# We define the function `T(p)` as

T(p) = TOCIM2

# (It turns out the circulation `T(p)` does not effectively depend on `p` but that's how we must define it anyway, i.e., as a function of `p`.)

# The local sources and sinks for the age take the form

function G(x,p)
    @unpack τ, z₀ = p
    return @. 1 - x / τ * (z ≤ z₀)
end

# as per the tracer equation.
# The `@unpack` line unpacks the parameters `τ` and `z₀`.
# The `return` line returns the net sources and sinks.
# (The `@.` "macro" tells Julia that the operations apply to every element.)

# We can define the vector `z` of depths with `depthvec`.

z = depthvec(grd)

# Now we must construct a type for `p` the parameters.
# This type must contain our parameters `τ` and `z₀`.

struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end

# The type is now ready for us to generate an instance of the parameter `p`.
# Let's use `τ = 1.0` (s) and `z₀` the minimum depth of the model.

p = IdealAgeParameters(1.0, 30.0)

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

plothorizontalslice(age_in_yrs, grd, depth=2000u"m", color=:magma)

# Or look at the horiontal mean

plothorizontalmean(age_in_yrs, grd)

# That's it for this tutorial...
# Good job!

