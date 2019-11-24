
#---------------------------------------------------------
# # Ideal mean age
#---------------------------------------------------------

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`ideal_mean_age.ipynb`](@__NBVIEWER_ROOT_URL__/examples/generated/ideal_mean_age.ipynb)

#---------------------------------------------------------
# ## The model
#---------------------------------------------------------

# We will simulate the ideal mean age of water.
# That is, the average amount of time since a water parcel had last contact with the surface.

# ### Tracer equation
#
# The ideal mean age is transported with water, is equal to $0$ at the surface, and increases by one second every second everywhere.
# In other words, the 3D field of the age, $a$, is governed by the tracer equation
#
# $$\frac{\partial a}{\partial t} + \nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right] a = 1,$$
#
# where $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right]$ is a differential operator that represents the transport by the ocean circulation.
# ($\boldsymbol{u}$ is the 3D vector field for the advection and $\mathbf{K}$ is the diffusivity matrix.)
# In the equation above, we also assume that there is the boundary condition that $a=0$ at the surface.

#---------------------------------
# ### Discretized tracer equation
#---------------------------------

# In AIBECS, the linear differential operator defined by $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right]$ is approximated by a constant matrix $\mathbf{T}$ when discretizing the continuous 3D ocean onto the model grid.
# This matrix can be small (e.g., for models with a few boxes), or large, like for the OCIM (more on the OCIM later).
# Similarly, the continuous 3D field of the age, $a$, is discretized into a column-vector, $\boldsymbol{a}$.
# (We represent scalars in italic, vectors in bold italic, and matrices in upstraight bold.)

# In the discrete case, we replace the boundary condition (that $a = 0$ at the surface) by imposing $\boldsymbol{a} = 0$ in the surface layer of the model grid.
# In practice, this is done by restoring $\boldsymbol{a}$ to $0$ with a very short timescale.
# The tracer equation thgus takes the form of
#
# $$\frac{\partial\boldsymbol{a}}{\partial t} = -\mathbf{T} \, \boldsymbol{a} + 1 - \boldsymbol{\Lambda} \, \boldsymbol{a},$$
#
# where $\boldsymbol{\Lambda}$ is a diagonal matrix with entries equal to $1 / \tau$ in the surface layer of the model grid and 0 otherwise.
# The timescale $\tau$ is chosen to be very small, ensuring that $\boldsymbol{a}$ is very close to $0$ at the surface.
# The first term represents the transport by the ocean circulation, the second term the source of 1 second per second everywhere, and the last term the fast relaxation.

#---------------------------------
# ### Steady-state
#---------------------------------

# The steady-state is the equilibrium that would be reached if we wait long enough for $a$ to not change anymore.
# Mathematically, the steady-state is also the state for which
#
# $$\frac{\partial a}{\partial t} = 0.$$

# Computationally, in the discrete case, this means that we just need to solve
#
# $$0 = -\mathbf{T} \, \boldsymbol{a} + 1 - \boldsymbol{\Lambda} \, \boldsymbol{a}$$
#
# to find $\boldsymbol{a}$.
# More specifically, we need to solve
#
# $$(\mathbf{T} + \boldsymbol{\Lambda}) \, \boldsymbol{a} = 1.$$

# Now that we have the equations laid down, let us chose the circulation transport matrix, $\mathbf{T}$.

#---------------------------------------------------------
# ## Using AIBECS
#---------------------------------------------------------

#md # !!! note
#md #     If this is the first time you are trying AIBECS, make sure you go through the prerequisites!
#nb # > **Note**
#nb # > If this is the first time you are trying AIBECS, make sure you go through the prerequisites!

# AIBECS can interpret tracer equations as long as you arrange them under the generic form:
#
# $$\frac{\partial \boldsymbol{x}}{\partial t} = \boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}),$$
#
# where $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p})$ is the rate of change of the state and $\boldsymbol{p}$ is the vector of model parameters.
# We only track the age here, so that the entire state of the system is determined by the age itself.
# In other words, here, $\boldsymbol{x} = \boldsymbol{a}$.

# We will use AIBECS to find the steady-state of the system.
# For AIBECS, this translates into finding the solution of $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0$.

# We start by telling Julia that we want to use the AIBECS package via

using AIBECS

#md # !!! note
#md #     If it's the first time you are running this line, the package will need precompiling.
#md #     This may take a minute or two.
#md #     (Just be patient... Or read on while you wait!)
#
#nb # > **Note**
#nb # > If it's the first time you are running this line, the package will need precompiling.
#nb # > This may take a minute or two.
#nb # > (Just be patient... Or read on while you wait!)
#
#md # !!! note
#md #     You should see a `Warning` for the `Flatten` package — just disregard it...
#md #     If you get an error though, please send me a copy of the output/error message, and I will try to troubleshoot it.
#nb # > **Note**
#nb # > You should see a `Warning` for the `Flatten` package — just disregard it...
#nb # > If you get an error though, please send me a copy of the output/error message, and I will try to troubleshoot it.



#---------------------------------
# ### The circulation
#---------------------------------

# We will use the circulation output from the Ocean Circulation Inverse Model (OCIM1).
# Basically, the OCIM provides researchers and oceanographers with a big sparse matrix that represents the global ocean circulation (advection and diffusion), which allows them to efficiently simulate the transport of passive tracers, like the age.
# (For more details, see Tim DeVries's [website](https://tdevries.eri.ucsb.edu/models-and-data-products/) and references therein.)
# With AIBECS, the OCIM0.1 and OCIM1 circulations can be loaded really easily, by simply typing

grd, T_OCIM = AIBECS.OCIM1.load()
typeof(T_OCIM), size(T_OCIM)

#
#md # !!! note
#md #     Julia may ask you to download the OCIM matrix for you, in which case you should say yes (i.e., type `y`).
#md #     Once downloaded, AIBECS will remember where it downloaded the file and it will only load it from your laptop.
#nb # > **Note**
#nb # > Julia may ask you to download the OCIM matrix for you, in which case you should say yes (i.e., type `y`).
#nb # > Once downloaded, AIBECS will remember where it downloaded the file and it will only load it from your laptop.
#
# Additionally to downloading the OCIM file, the `load()` command loads 2 variables in the Julia workspace:
# - `grd` — a `OceanGrid` object containing information about the 3D grid of the OCIM circulation, like the latitude, longitude, and depth of each grid boxes, as well as information on which boxes are "wet" or "dry".
# - `T_OCIM` — the transport matrix representing advection and diffusion.
#
# The second line in command above tells you the type and the size of `T_OCIM`.
# It is a sparse matrix (CSC just means that it is stored in Compressed Sparse Column format) and is quite big!
#md # !!! note
#md #     A sparse matrix is just a matrix with very few non-zero entries.
#md #     Computationally, sparse matrices are stored differently than full matrices to save memory (no need to save all those zeros), and are much faster to use too!
#nb # > **Note**
#nb # > A sparse matrix is just a matrix with very few non-zero entries.
#nb # > Computationally, sparse matrices are stored differently than full matrices to save memory (no need to save all those zeros), and are much faster to use too!

# Anyway, this looks good, so let's move on with setting up the model!
#
# We have already loaded the transport matrix, `T_OCIM`, for the ocean circulation, but we must tell AIBECS that it applies to the age.
# To do that, we define a function of the parameters (although there are no parameters involved in this case, this is just the way AIBECS works for the moment).

T_age(p) = T_OCIM

# (Functions in Julia can be created in one line, just as above.)
# That's it for the circulation.
# Now, let's define the local sources and sinks.

#---------------------------------
# ### The local sources and sinks
#---------------------------------

# We will denote the age, $\boldsymbol{a}$, by the variable `age` in Julia.
# (It's good practice to use explicit names!)
# We need to translate the local sources and sinks in our discretized state function $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p})$ into Julia code.

# #### The source
#--------------------------

# Remember the age increases by $1$ second every second and everywhere.
# So its source function is equal to, well, `1`! (seconds per seconds means it is unitless).
# Let's create the local source function:

source_age(age, p) = 1

# #### The sink
#--------------------------

# Recall that the age must also be $0$ at the surface.
# And that we implement this boundary condition by restoring the age very fast to $0$ in the surface layer.
# This will act as the sink for the age.
# But first, we need to figure out where "the surface layer" is.
# To do that, AIBECS can generate a number of useful constants for you.
# (You can see the list of functions by typing `varinfo(AIBECS)` at the REPL.)
# Here we will use the vector of grid box depths, `z`, which AIBECS can generate for us via

z = vector_of_depths(grd)

# So what is the top layer?
# Let's investigate what's the minimum depth:

minimum(z)

# The surface layer in the OCIM grid has its center at about $18\,$m depth.
# We can create a mask of the surface layer via `z .< 20`.
# (This will return a vector of `0`s and `1`s, depending on whether the depth, `z`, is less than `20`.)
#
#md # !!! note
#md #     In Julia (like in MATLAB), placing a dot, `.`, in front of operators is a convenient way to do element-wise operations.
#nb # > **Note**
#nb # > In Julia (like in MATLAB), placing a dot, `.`, in front of operators is a convenient way to do element-wise operations.
#
# Then, we implement the local sink by restoring the age to `0` with a timescale `τ`, via

function sink_age(age, p)
    τ = p.τ
    return age .* (z .< 20u"m") / τ
end

#md # !!! note
#md #     Julia allows you to use unicode for your functions and variables, like for `τ`.
#nb # > **Note**
#nb # > Julia allows you to use unicode for your functions and variables, like for `τ`.
#
# Here, we have defined a Julia function using the `function` keyword because the sink is a bit more complicated, so that we needed two lines to define it.
# The first line unpacks the model parameters, which is just the restoring timescale, `τ`, in this case.
# We will chose the value for `τ` later.

# #### Net sources and sinks
#--------------------------

# The sources minus the sinks are simply defined by

sms_age(age, p) = source_age(age, p) .- sink_age(age, p)

# #### Model parameters
#--------------------------

# We must define the parameters... And AIBECS comes with an API for that!

t = empty_parameter_table()    # initialize table of parameters
add_parameter!(t, :τ, 1u"s")   # add the parameter we want (τ = 1s)
initialize_Parameters_type(t, "IdealAgeParameters")  # Generate the parameter type
t

# Note, in particular, that we gave our parameter `τ` a unit.
# Yes, Julia comes with some nice functionality to deal with units directly!
# The lines above created a table that contains all the info for generating the parameters vector, $\boldsymbol{p}$.
# To generate the parameters in AIBECS we do:

p = IdealAgeParameters()

# where we have used the constructor `IdealAgeParameters`, whose name we defined in the previous cell.
# Here we did not really need to create `p` as a parameters vector, since it has only one element, `τ`, in it.
# However, we are here to learn, and this structure and functionality comes in very handy when one deals with many parameters.
# (And as you can imagine, having all the parameters in a nice table ready for being used in a publication comes quite handy!)



# #### State function and Jacobian
#--------------------------

# Similarly to `p`, let's create a state `x` to start with.
# The vector `x` will be our initial guess for the state.
# Let's assume that the age is `1` (seconds) everywhere (as an initial guess):

nb = number_of_wet_boxes(grd)  # number of wet boxes
x = ones(nb)

# The first line above defines the number of wet grid boxes, `nb`.
# Here, this is also the length of the state vector `x`, because there is only one tracer, `age`.
# In the second line, the `ones` function creates a vector of `1`s of the size you give it (the number of wet grid boxes, `nb`, here, which we defined earlier).

# Finally, the last step for the set up is to define $\boldsymbol{F}$.
# Using AIBECS, this is done via

T_matrices = (T_age,)           # bundles all the transport matrices in a tuple
sources_minus_sinks = (sms_age,) # bundles all the source-sink functions in a tuple
F, ∇ₓF = state_function_and_Jacobian(T_matrices, sources_minus_sinks, nb) # generates the state function (and its Jacobian!)
F(x,p)

# That's it!
# We have just created a model of the mean age.
#
#
# Lines 2 and 3 are just telling AIBECS
# - what transport matrices it should use for the transport of these tracers, and
# - and what local sources and sinks should be appplied to these tracers
#
#md # !!! note
#md #     The `(x,)` syntax returns a tuple of one element — the comma is necessary because without it, `(x)` would be just like `x` with brackets around it.
#md #     This interface of AIBECS was developed for case with multiple tracers in mind, and might look a bit odd for a single tracer.
#md #     But in the future, this might be cleaned up to be easier to work with single tracers.
#nb # > **Note**
#nb # > The `(x,)` syntax returns a tuple of one element — the comma is necessary because without it, `(x)` would be just like `x` with brackets around it.
#nb # > This interface of AIBECS was developed for case with multiple tracers in mind, and might look a bit odd for a single tracer.
#nb # > But in the future, this might be cleaned up to be easier to work with single tracers.

# The fourth line creates two functions:
# - `F` — the numerical version of the **state function**, $\boldsymbol{F}$, of our model of the mean age, and
# - `∇ₓF` — the **Jacobian matrix** of the state function, i.e., $\nabla_{\boldsymbol{x}}\boldsymbol{F}$.
# Yes, AIBECS just automatically created an exact derivative of your input, using autodifferentiation via dual numbers.
# (I'd be very excited to detail how this is implemented here, but it is an entirely different discussion.)

# The last line just checks that our generated `F` works with our initial guess `x` and parameter vector `p`.

# #### Solving for the steady-state
#--------------------------

# The Jacobian, `∇ₓF` is essential to solving the steady-state equation $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0$ fast.
# Specifically, solving $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0$ is done via Newton's method.
# By starting from an initial guess, that you will have to provide, it will iterate over this recursion relation
#
# $$\boldsymbol{x}_{k+1} = \boldsymbol{x}_{k} - \nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x}_{k}, \boldsymbol{p})^{-1} \boldsymbol{F}(\boldsymbol{x}_{k}, \boldsymbol{p})$$
#
# until $\boldsymbol{F}(\boldsymbol{x}_{k}, \boldsymbol{p})$ is sufficiently small.
# Now I should note that here, our age model is linear in $x$ (or `age` in our code), so that the solution will be found in a single iteration, or a sinle "matrix inversion", as could be seen from our steady-state equation for $\boldsymbol{a}$.

#md # !!! note
#md #     AIBECS comes with a built-in algorithm and an API to solve for the steady-state, so you don't have to worry about all these details!
#nb # > **Note**
#nb # > AIBECS comes with a built-in algorithm and an API to solve for the steady-state, so you don't have to worry about all these details!



# ##### Define the Steady-state problem in AIBECS

# First, we create an instance of the steady-state problem, via

prob = SteadyStateProblem(F, ∇ₓF, x, p)

# where we have simply provided the state function, $\boldsymbol{F}$, the Jacobian, $\nabla_{\boldsymbol{x}}\boldsymbol{F}$, the initial guess and the parameters.
# The `SteadyStateProblem` function is a standard "DiffEqBase" constructor that I have overloaded in my package so that you can easily generate the model here.



# ##### Solve for the steady-state with AIBECS

# Finally, we can find the solution in litterally one line, via the `solve` function:

age = solve(prob, CTKAlg())

# Here, I have provided the `solve` function with two things:
# - the problem, `prob`, which we just defined, and
# - the quasi-Newton algorithm that I wrote in Julia, denoted by `CTKAlg()` after C.T. Kelley, who originally wrote it in MATLAB.
#
# The last line should take about 10 seconds to 1 minute, depending on your laptop.
# That's it!
# We solved for the steady state!
# Everyone here deserves a nice tap on the shoulder — Good job!
# Now let's see what this age looks like on a map



#---------------------------------------------------------
# ## Figures
#---------------------------------------------------------

# First, let's tell Julia we want to use Plots.jl.

using Plots

# We will plot the age at a depth of about 1000m.
# For convenience, we would rather express the age in years
# (instead of the default SI unit of time, seconds, which is too small).

age_in_yrs = age * u"s" .|> u"yr"

# A few things happened here.
# - `age * u"s"` attaches the seconds unit to it
# - `.|> u"yr"` converts it to years (the dot in front ensures that the conversion is done for every box

horizontalslice(age_in_yrs, grd, 2000; color=:magma)

# Note how the plotting function `horizontalslice` knew about the unit of the age and placed it on the colorbar! Nice, right?

#md # !!! tip
#md #     The function `horizontalslice(x, grd, depth)` is provided by AIBECS.
#md #     It is actually just a recipe for Plots.jl that figures out a few things for you,
#md #     like extracting the slice at the given depth.
#md #     But you can customize it to your liking, by appending keyword arguments,
#md #     like `color=:magma` here.
#md #     Head over to the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package
#md #     documentation to see a more complete list of attributes

#nb # > **Tip!**
#nb # > The function `horizontalslice(x, grd, depth)` is provided by AIBECS.
#nb # > It is actually just a recipe for Plots.jl that figures out a few things for you,
#nb # > like extracting the slice at the given depth.
#nb # > But you can customize it to your liking, by appending keyword arguments,
#nb # > like `color=:magma` here.
#nb # > Head over to the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package
#nb # > documentation to see a more complete list of attributes


# At 1000m, the age ranges from a few years below deep water formation regions (Wedell Sea, North Atlantic), and reaches a dozen of centuries in the North Pacific!
# Anyway, that's pretty much it!
# Good job!

