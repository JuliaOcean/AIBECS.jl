# ## Ideal mean age
#
# **Prerequisites**: If this is the first time you are using AIBECS, make sure you have installed the required packages required this example:
# - [AIBECS](https://github.com/briochemc/AIBECS.jl)
# - [PyPlot.jl (update link)]()
# - [PyCall.jl (update link)]()
# - [Conda.jl (update link)]() + `conda install cartopy`
#
# We start by making sure Julia knows that we are going to use the AIBECS package:

using AIBECS

# If it's the first time you are running this line, the package will need precompiling.
# This may take a minute or two.
# Just be patient :)
# You should see a `Warning` for the `Flatten` package — just disregard it...
# If you get an error though, please send me a copy of the output/error message, and I will try to troubleshoot it.

# ### Set up the model
#
# Let's start with a very simple model.
# We will simulate the ideal mean age of water.
# That is, the average amount of time since a water parcel had last contact with the surface.
# Here the "surface" will be the top layer of the model grid.
# So, we have one tracer, the age, which we will denote by `age` in Julia.
# But before diving into the equations, let us chose the circulation.

# #### The circulation
#
# We will use the circulation output from the Ocean Circulation Inverse Model (OCIM) version 1.1.
# Basically, the OCIM provides scientists, like you and me, with a big sparse matrix that represents the global ocean circulation (advection and diffusion).
# (For more details, see Tim DeVries's [website](https://tdevries.eri.ucsb.edu/models-and-data-products/) and references therein.)
# With AIBECS, the OCIM circulation can be loaded really easily, by simply typing

const wet3d, grd, T_OCIM = AIBECS.OCIM1.load() ;

# Here, I have used the `const` prefix to tell Julia these are constants (it helps Julia's compiler)
# I also use a semicolon, `;`, to suppress the output.
# Anyway, Julia may ask you to download the OCIM matrix for you, in which case you should say yes (i.e., type `y`).
# Once downloaded, AIBECS will remember where it downloaded the file and it will only load it from your laptop.
# Additionally to downloading the OCIM file, the last command loads 3 variables in the Julia workspace:
# - `grd` contains information about the 3D grid of the OCIM circulation, like the latitude, longitude, and depth of each grid boxes. `grd` is a dictionary (a `Dict` in Julia, equivalent to a `struct` in MATLAB).
# - `wet3d` is a 3D array of the model grid, filled with `1`'s at "wet" grid boxes and `0`'s and "land" grid boxes.
# - `T_OCIM` is the transport matrix representing advection and diffusion.
#
# Let's check that `T_OCIM` is what we think it is, by using the `typeof` command (which tells you, well, the "type" of things):

typeof(T_OCIM)

# Great, it's a sparse matrix! (CSC just means that it is stored in Compressed Sparse Column format.)

# We can get its size with the `size` command:

size(T_OCIM)

# Yes, it's quite a big matrix!
# (But it is sparse, so it has very few non-zero entries.)
#
# Anyway, this looks good, so let's move on with setting up the model!
#
# The state of the model, $\boldsymbol{x}$, is entirely defined by one tracer: the age (denoted by the `age` variable).
# The age is transported along with water parcels, so its transport matrix is the ocean circulation matrix that we just loaded: `T_OCIM`.
# We tell AIBECS that by first writing the transport matrix for the age as a function of the parameters (although there are no parameters in this simple model).

T_age(p) = T_OCIM

# (Functions in Julia can be created in one line, just as above.)

# #### The local sources and sinks
#
# Now what are the local sources and sinks of `age`?
# The age increases by $1$ second every second and everywhere.
# So its source is $1$ (seconds per seconds means it is unitless) everywhere.
# Let's create the local source function:

source_age(age, p) = 1

# The age is also $0$ at the surface.
# To enforce that, we restore the age very fast to $0$ at the surface.
# This will act as the sink for the age.
# But first, we need to figure out where "the surface" is.
# To do that, AIBECS can generate a number of useful constants for you, like the vector of depths or the vector of volumes.
# These are computed this way

const iwet = indices_of_wet_boxes(wet3d)
const nb = number_of_wet_boxes(wet3d)
const v = vector_of_volumes(wet3d, grd)
const z = vector_of_depths(wet3d, grd)
const ztop = vector_of_top_depths(wet3d, grd)

# Here, we will only use `iwet`, `nb`, and `z`.
# Just for your informtation, these are:
# - `iwet` — the vector of the indices of the wet boxes, i.e., those boxes from the 3D ocean that are, well, not land.
# - `nb` — the number of wet grid boxes. Here, this is also the length of the state vector `x`, because there is only one tracer, `age`.
# - `v` — the vector of grid box volumnes.
# - `z` — the vector of grid box depths.
# - `ztop` — the vector of the depths of the top of the grid boxes.

# So what is the top layer?
# Let's investigate what's the minimum depth:

minimum(z)

# So the surface layer in the OCIM grid has its center at about 18m depth.
# We can create a mask of the surface layer via `z .< 20`.
# (This will return a vector of `0`s and `1`s, depending on whether the depth, `z`, is less than `20`.)
#
# Then, we implement the local sink by restoring the age to `0` with a timescale `τ`, via

function sink_age(age, p)
    τ = p.τ
    return age .* (z .< 20) / τ
end

# (Julia allows you to use unicode for your functions and variables, like for `τ`.)
# Here, we have defined a Julia function using the `function` keyword because the sink is a bit more complicated, so that we needed two lines to define it.
# The first line unpacks the model parameters, which is just `τ` in this case.
# We will chose the value for `τ` later.
#
# In the function `sink_age`, you might notice that we used a `.*` and `.<` instead of simply `*` and `<`.
# This is called "broadcast" in Julia, and it means it is an element-wise operation.
# Hence, `sink_age` will be large wherever `age` is large, but only at the surface, where `z .< 20`, which is exactly what we want.
#
# Now the sources minus the sinks is simply created by

sms_age(age, p) = source_age(age, p) .- sink_age(age, p)

# where we have used `.-` instead just `-` (to "broadcase" the age source, which is just the scalar `1`).
#
# Finally, the last step for the set up is to define $\boldsymbol{F}$.
# Using AIBECS, this is done via

T_matrices = (T_age,)           # bundles all the transport matrices in a tuple
sources_minus_sinks = (sms_age,) # bundles all the source-sink functions in a tuple
F, ∇ₓF = state_function_and_Jacobian(T_matrices, sources_minus_sinks, nb) # generates the state function (and its Jacobian!)

# That's it!
# We have just created a model of the mean age.
# The first 2 lines in the cell above are just telling AIBECS
# - what transport matrices it should use for the transport of these tracers, and
# - and what local sources and sinks should be appplied to these tracers
#
# > Side note: This interface was developed for multiple tracers, and might look a bit odd for a single tracer.
# > (The `(x,)` syntax returns a tuple of one element — the comma is necessary because without it, `(x)` would be just like `x` with brackets around it.)

# The last line creates two functions:
# - `F` — the numerical version of the **state function**, $\boldsymbol{F}$, of our model of the mean age, and
# - `∇ₓF` — the **Jacobian matrix** of the state function, i.e., $\nabla_{\boldsymbol{x}}\boldsymbol{F}$.
# Yes, AIBECS just automatically created an exact derivative of your input, using autodifferentiation via dual numbers.
# (I'd be very excited to detail how this is implemented here, but it is an entirely different discussion.)
#
# The Jacobian, `∇ₓF` is essential to solving the steady-state equation $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0$ fast.
# Specifically, solving $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0$ is done via Newton's method.
# By starting from an initial guess, that you will have to provide, it will iterate over this recursion relation
#
# $$\boldsymbol{x}_{k+1} = \boldsymbol{x}_{k} - \nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p})^{-1} \boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) $$
#
# until if finds a solution.
# Now I should note that here, our age model is linear in $x$ (or `age` in our code), so that the solution will be found in a single iteration, or a sinle "matrix inversion".

# Anyway, AIBECS is nice, and comes bundled with an algorithm and an API to solve for the steady-state, so you don't have to worry about all these cumbersome details!
# But before we do so, we must first define the parameter... And AIBECS also comes with an API for that!

t = empty_parameter_table()    # initialize table of parameters
add_parameter!(t, :τ, 1u"s")   # add the parameter we want (τ = 1s)
initialize_Parameters_type(t)  # Generate the parameter type

# The lines above created a table that contains all the info for my parameters, and a corresponding parameter vector, `p₀`.
# Note, in particular, that we gave our parameter `τ` a unit.
# Yes, Julia comes with some nice functionality to deal with units directly!
#
# Let us have a look at the table `t` and the parameter `p₀`

t

p₀

# Here we did not really need to create `p₀` as a parameters vector, since it has only one element, `τ`, in it.
# However, we are here to learn, and this structure and functionality comes in very handy when one deals with many parameters.
# (And as you can imagine, having all the parameters in a nice table ready for being used in a publication comes quite handy!)
#
# Anyway, now let's finish our model, and set an initial guess for the age.
# Let's assume the age is `1` (seconds) everywhere ():

x₀ = ones(nb)

# the `ones` function creates a vector of `1`s of the size you give it (the number of wet grid boxes, `nb`, here, which we defined as a constant earlier).

# ### Run the simulation, i.e., solve for the steady state
#
# First, we create an instance of the steady-state problem, via

prob = SteadyStateProblem(F, ∇ₓF, x₀, p₀)

# where we have simply provided the state function, $\boldsymbol{F}$, the Jacobian, $\nabla_{\boldsymbol{x}}\boldsymbol{F}$, the initial guess and the parameters.
# The `SteadyStateProblem` function is a standard "DiffEqBase" constructor that I have overloaded in my package so that you can easily generate the model here.
#
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

# ### Make some figures
#
# We will plot a horizontal slice of the age at about 1000m depth.

# First, we must take `age`, and rearrange it into the 3D grid.
# For that we use `iwet`, which we defined as a constant above (recall `iwet` is the vector of the indices of wet points in the 3D grid).
#
# We can now rearrange `age` in these few lines of code:

age_3D = NaN * wet3d # creates a 3D array of NaNs
age_3D[iwet] = age   # Fills the wet grid boxes with the age values
size(age_3D)         # Just to check the size of age_3D

# Thus, `age_3D` is a 3D-array of which I just showed you the size, which corresponds to the OCIM grid.
#
# Now we must find the index of the depth that is closest to 1000m.
# To do that we must use the depth information contained in `grd`.
# Let us first create a small vector of the depths of the grid:

depth = vec(grd["zt"])

# So we could count the index of the entry we want, or we could use the `findfirst` function, like below.
# (Feel free to change the value of `iz` if you want to see a slice at another depth.)

iz = findfirst(depth .> 1000)
iz, depth[iz]

# So the 13-th layer lies at 1104m, which should do the trick for us.
#
# Finally, we need the latitude and longitudes of the grid, contained in `grd`.
# As for `depth`, we can use the OCIM's `grd` output:

lat, lon = vec(grd["yt"]), vec(grd["xt"])

# So these are the latitudes and longitudes of the map we are about to plot.

# A last thing we can do is convert the age from seconds, `u"s"`, to years, `u"yr"`, thanks to the Unitful package 
# (loaded automatically by AIBECS).

age_3d_1000m_yr = age_3D[:,:,iz] * ustrip(1.0u"s" |> u"yr")

# Finally! Let's have a look at this ideal mean age!
# To make figures, we will use Cartopy (that you should have `add`ed if you went through the README correctly).
# To use it we simply type

ENV["MPLBACKEND"]="qt5agg"
using PyPlot, PyCall

# The first line is needed for Mac users. It's a bug that should eventually be resolved, but for now this seems to make it work.
# 
# We then import cartopy, define a new plot, add some coastlines because they are pretty, and add our slice of age at 1000m depth to it via

ccrs = pyimport("cartopy.crs")
ax = subplot(projection=ccrs.Robinson(central_longitude=-155.0))
ax.coastlines()
# making it cyclic for Cartopy
lon_cyc = [lon; 360+lon[1]] 
age_cyc = hcat(age_3d_1000m_yr, age_3d_1000m_yr[:,1])
# And plot
p = contourf(lon_cyc, lat, age_cyc, levels=0:100:1200, transform=ccrs.PlateCarree(), zorder=-1)
colorbar(p, orientation="horizontal")

# That's it!
# Good job!
# 
# At 1000m, the age ranges from a few years below deep water formation regions (Wedell Sea, North Atlantic), and reaches a dozen of centuries in the North Pacific! 
# This is pretty good for so little work!