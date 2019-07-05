
#---------------------------------------------------------
# # Optimization of a simple phosphorus cycling model
#---------------------------------------------------------

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`P_model_2_tracers.ipynb`](@__NBVIEWER_ROOT_URL__examples/generated/P_model_2_tracers.ipynb)

#md # !!! warn
#md #     This example is not finished yet.
#md #     It lacks a number of things, including explanatory text, references, and more.
#md #     However, it can already serve as an example of model optimization.

#nb # > **Warning!**
#nb # > This example is not finished yet.
#nb # > It lacks a number of things, including explanatory text, references, and more.
#nb # > However, it can already serve as an example of model optimization.

#---------------------------------------------------------
# ## Tracer equations
#---------------------------------------------------------

# We consider a simple model for the cycling of phosphorus with 2 state variables consisting of phosphate (PO₄) AKA dissolved inorganic phosphorus (DIP) and particulate organic phosphorus (POP).
# The dissolved phases are transported by advection and diffusion whereas the particulate phase sinks rapidly down the water column without any appreciable transport by the circulation.
#
# The governing equations that couple the 3D concentration fields of DIP and POP, denoted $x_\mathsf{DIP}$ and $x_\mathsf{POP}$, respectively, are:
#
# $$\left[\frac{\partial}{\partial t} + \nabla \cdot (\boldsymbol{u} + \mathbf{K}\cdot\nabla )\right] x_\mathsf{DIP} = -U(x_\mathsf{DIP}) + R(x_\mathsf{POP}),$$
#
# and
#
# $$\left[\frac{\partial}{\partial t} + \nabla \cdot \boldsymbol{w}\right] x_\mathsf{POP} = U(x_\mathsf{DIP}) - R(x_\mathsf{POP}).$$

# The 3D vector $\boldsymbol{u}$ is the fluid velocity and $\mathbf{K}$ is the eddy-diffusion tensor.
# Thus, $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right]$ is a differential operator that represents the transport by the ocean circulation.
# The 3D vector $\boldsymbol{w}$ is the particulate sinking velocity, so that the operator $\nabla \cdot \boldsymbol{w}$ acting on $x_\mathsf{POP}$ is the divergence of the flux of POP.
# We will assume that it increases linearly with depth, i.e., that
#
# $$\boldsymbol{w}(\boldsymbol{r}) = \boldsymbol{w}_0 + \boldsymbol{w}' z,$$
#
# where $\boldsymbol{w}_0$ and $\boldsymbol{w}'$ are simple scalars, independent of location $\boldsymbol{r}$.
# ($z$ represents the depth of location $\boldsymbol{r}=(x,y,z)$.)

# The function $U$ represents the biological uptake of DIP by phytoplankton, which we model as
#
# $$U(x_\mathsf{DIP}) = \frac{x_\mathsf{DIP}}{\tau} \, \frac{x_\mathsf{DIP}}{x_\mathsf{DIP} + k} \, (z < z_0),$$
#
# where the timescale, $\tau$, sets the maximum specific DIP uptake rate.
# The first term in the product scales linearly with DIP, which is a nice property to have in general.
# It also makes it easier to reason about the extreme cases.
# However, here the uptake is not linear:
# the second term in the product is a hyperbolic term (after Monod or Michaelis-Menten) and reduces the slope of the uptake function when DIP is low.
# When DIP is very low, the slope is close to 0.
# But when DIP is very high, the slope is close to $1/τ$.
# The last term is to be understood as a condition, i.e., that the uptake happens only in the euphotic layer, where there is sufficient light for photosynthesis to occur.
# The depth $z_0$ is the depth of the euphotic layer.

# The function $R$ defines the remineralization rate of POP, which converts POP back into DIP.
# For the remineralization, we simply use a linear rate constant, i.e.,
#
# $$R(x_\mathsf{POP}) = \kappa \, x_\mathsf{POP}.$$



#---------------------------------------------------------
# ## Translation into AIBECS code
#---------------------------------------------------------

# ### Transport matrices

# We start by telling Julia we want to use the AIBECS

using AIBECS

# We load the circulation and grid of the OCIM0.1 via

wet3D, grd, T_OCIM = OCIM0.load() ;

# where the operator $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right]$ is represented by the sparse matrix `T_OCIM`.
# In AIBECS, all the transport matrices must be declared as functions of the parameters `p`, so we simply define

T_DIP(p) = T_OCIM

# For the sinking of particles, we are going to create a sparse matrix that depends on the parameters, too, which define how fast particles sink.

T_POP(p) = buildPFD(grd, wet3D, sinking_speed = w(p))

# for which we need to define the sinking speed `w(p)` as a function of the parameters `p`.
# Following the assumption that it increases linearly with depth, we write it as

w(p) = p.w₀ .+ p.w′ * z

# For this to work, we must create a vector of depths, `z`, which is simply done via

iwet = findall(wet3D)
z = ustrip.(grd.depth_3D[iwet])

# (We must strip the units via `ustrip` because units do not percolate through the functionality of AIBECS flawlessly at this early stage of development.)
# We combine the transport operators for DIP and POP into a single tuple via

T_all = (T_DIP, T_POP) ;



# ### Local sources and sinks

# ##### Geological Restoring

# Because AIBECS will solve for the steady state solution directly without time-stepping the goverining equations to equilibrium, we don't have any opportunity to specify any intial conditions.
# Initial conditions are how the total amount of conserved elements get specified in most global biogeochemical modelels.
# Thus to specify the total inventory of P in AIBECS we add a very weak resporing term to, e.g., the DIP equation.
# The time-scale for this restoring term is chosen to be very long compared to the timescale with which the ocean circulation homogenizes a tracer.
# Because of this long timescale we call it the geological restoring term, but geochemists who work on geological processes don't like that name!
# In any event the long timescale allows us to prescribe the total inventory of P in a way that yields the same solution we would have gotten had we time-stepped the model to steady-state with the total inventory prescribed by the initial condition.
# We define the geological restoring simply as

geores(x, p) = (p.xgeo .- x) / p.τg

# ##### Uptake (DIP → POP)

# For the phosphate uptake, we must take some precautions to avoid the pole of the hyperbolic term.
# That is, we ensure that uptake only occurs when DIP is positive via the `relu` function

relu(x) = (x .≥ 0) .* x

# Then we simply write the uptake as a Julia function, following the definition of $U$:

function uptake(DIP, p)
    τ, k, z₀ = p.τ, p.k, p.z₀
    DIP⁺ = relu(DIP)
    return 1/τ * DIP⁺.^2 ./ (DIP⁺ .+ k) .* (z .≤ z₀)
end

# where we have "unpacked" the parameters to make the code clearer and as close to the mathematical equation as possible.

# ##### Remineralization (POP → DIP)

remineralization(POP, p) = p.κ * POP

# ##### Net sources and sinks

# We lump the sources and sinks into `sms` functions (for Sources Minus Sinks):

function sms_DIP(DIP, POP, p)
    return -uptake(DIP, p) + remineralization(POP, p) + geores(DIP, p)
end
function sms_POP(DIP, POP, p)
    return  uptake(DIP, p) - remineralization(POP, p)
end

# and create a tuple of those just like for the transport matrices via

sms_all = (sms_DIP, sms_POP) # bundles all the source-sink functions in a tuple


# ### Parameters

# We now define and build the parameters.
# In AIBECS, this is done by first creating a table of parameters, using the API designed for it.
# First, we create an empty table:

t = empty_parameter_table()    # empty table of parameters

# Then, we add every parameter that we need

add_parameter!(t, :xgeo, 2.12u"mmol/m^3",
               optimizable = true,
               variance_obs = ustrip(upreferred(0.1 * 2.17u"mmol/m^3"))^2,
               description = "Mean PO₄ concentration")
add_parameter!(t, :τg, 1.0u"Myr",
               description = "Geological restoring timescale")
add_parameter!(t, :k, 6.62u"μmol/m^3",
               optimizable = true,
               description = "Half-saturation constant (Michaelis-Menten)")
add_parameter!(t, :z₀, 80.0u"m",
               description = "Depth of the euphotic layer base")
add_parameter!(t, :w₀, 0.64u"m/d",
               optimizable = true,
               description = "Sinking velocity at surface")
add_parameter!(t, :w′, 0.13u"1/d",
               optimizable = true,
               description = "Vertical gradient of sinking velocity")
add_parameter!(t, :κ, 0.19u"1/d",
               optimizable = true,
               description = "Remineralization rate constant (POP to DIP)")
add_parameter!(t, :τ, 236.52u"d",
               optimizable = true,
               description = "Maximum uptake rate timescale")

# These values of these parameters have been optimized for the same model but embedded in the OCIM1 circulation and should yield a pretty good fit to observations.
# However, because the circulation has changed from the OCIM1 to the OCIM0.1 used here, the best parameters have probably changed, too.
# Thus, we will re-optimize these parameters below.
# But first, we will just try to find the steady-state solution with these initial parameter values.

# We generate a new type for the parameters, via

initialize_Parameters_type(t, "Pcycle_Parameters")   # Generate the parameter type

# Finally, we can instantiate the parameter vector in a single line:

p = Pcycle_Parameters()

# where we have used the name of the new type that we created just earlier.

#---------------------------------------------------------
# ## Finding the steady-state solution
#---------------------------------------------------------

# ### Generate state function and Jacobian

# We generate the state function and its Jacobian in a single line

nb = length(iwet)
F, ∇ₓF = state_function_and_Jacobian(T_all, sms_all, nb)

# ### Solve for the steady-state

# We start from an initial guess,

x = p.xgeo * ones(2nb) # initial iterate

# Create an instance of the steady-state problem

prob = SteadyStateProblem(F, ∇ₓF, x, p)

# and solve it

s = solve(prob, CTKAlg()).u

# ### Plotting

# For plotting, we first unpack the state

DIP, POP = state_to_tracers(s, nb, 2)

# We will plot the concentration of DIP at a given depth horizon

iz = findfirst(grd.depth .> 2000u"m")
iz, grd.depth[iz]

#-

DIP_3D = rearrange_into_3Darray(DIP, wet3D)
DIP_2D = DIP_3D[:,:,iz] * ustrip(1.0u"mol/m^3" |> u"mmol/m^3")
lat, lon = ustrip.(grd.lat), ustrip.(grd.lon)

# Create the Cartopy canvas

ENV["MPLBACKEND"]="qt5agg"
using PyPlot, PyCall
clf()
ccrs = pyimport("cartopy.crs")
ax = subplot(projection = ccrs.EqualEarth(central_longitude=-155.0))
ax.coastlines()

# Making the data cyclic in order for the plot to look good in Cartopy

lon_cyc = [lon; 360+lon[1]]
DIP_2D_cyc = hcat(DIP_2D, DIP_2D[:,1])

# And plot

plt = contourf(lon_cyc, lat, DIP_2D_cyc, levels=0:0.2:3.6, transform=ccrs.PlateCarree(), zorder=-1)
colorbar(plt, orientation="horizontal");
title("PO₄ at $(string(round(typeof(1u"m"),grd.depth[iz]))) depth using the OCIM0.1 circulation")
gcf()



#---------------------------------------------------------
# ## Optimizing the model parameters
#---------------------------------------------------------

# To optimize the model, we will use data from the World Ocean Atlas (WOA).
# Specifically, we will use the observed PO<sub>4</sub> climatology, which the WOA produces every few years.
# Thanks to the [WorldOceanAtlasTools](https://github.com/briochemc/WorldOceanAtlasTools.jl) package combined with some AIBECS functionality, this is done a just a few commands:

using WorldOceanAtlasTools
μDIPobs3D, σ²DIPobs3D = WorldOceanAtlasTools.fit_to_grid(grd, 2018, "phosphate", "annual", "1°", "an")
μDIPobs, σ²DIPobs = μDIPobs3D[iwet], σ²DIPobs3D[iwet]
μx = (μDIPobs, missing)
σ²x = (σ²DIPobs, missing)

# We define some hyper parameters for the weight of the mismatch of DIP, POP, and the parameters with observations / their prior estimates.

ωs = [1.0, 0.0] # the weight for the mismatch (weight of POP = 0)
ωp = 1e-4       # the weight for the parameters prior estimates

# We then generate the objective function and some derivatives using the PO₄ mean and variance

v = ustrip.(vector_of_volumes(wet3D, grd))
f   =   generate_objective(ωs, μx, σ²x, v, ωp, mean_obs(p), variance_obs(p))
∇ₓf = generate_∇ₓobjective(ωs, μx, σ²x, v, ωp, mean_obs(p), variance_obs(p))
∇ₚf = generate_∇ₚobjective(ωs, μx, σ²x, v, ωp, mean_obs(p), variance_obs(p))

# For the optimization, we will use the F-1 algorithm from the [F1Method](https://github.com/briochemc/F1Method.jl) package.

using F1Method

# For the F-1 algorithm, we must initialize the memory buffer

mem = F1Method.initialize_mem(s, p)

# We then define the objective, gradient, and Hessian using the F1Method package API

objective(p) = F1Method.objective(f, F, ∇ₓF, mem, p, CTKAlg())
gradient(p) = F1Method.gradient(f, F, ∇ₓf, ∇ₓF, mem, p, CTKAlg())
hessian(p) = F1Method.hessian(f, F, ∇ₓf, ∇ₓF, mem, p, CTKAlg())

# Because we want the parameters to remain positive, we define the change of variables
# The change variable from vector `λ = log(p)` to Parameters `p`:

λ2p(λ) = AIBECS.opt_para(exp.(λ))

# and the corresponding derivatives:

∇λ2p(λ) = exp.(λ)
∇²λ2p(λ) = exp.(λ)

# Reverse change of variables from Parameters p to vector λ

p2λ(p) = log.(optvec(p))
λ = p2λ(p)

obj(λ) = objective(λ2p(λ))
using LinearAlgebra
grad(λ) = gradient(λ2p(λ)) * Diagonal(∇λ2p(λ))
function hess(λ)
    ∇p = Diagonal(∇λ2p(λ)) # for variable change
    ∇²p = Diagonal(∇²λ2p(λ)) # for variable change
    G = vec(gradient(λ2p(λ)))
    H = hessian(λ2p(λ))
    return ∇p * H * ∇p + Diagonal(G) * ∇²p
end

# Use optim to optimize

using Optim

# Add storage for use with Optim

m = length(p)
grad(s, λ) = s[1:m] .= vec(grad(λ))
hess(s, λ) = s[1:m,1:m] .= hess(λ)

# Set options of Optim

opt = Optim.Options(store_trace = false, show_trace = true, extended_trace = false)

# Run optimization

results = optimize(obj, grad, hess, λ, NewtonTrustRegion(), opt)

# Because we started from an other optimized set of parameters, this should only take a few iterations.
# The starting parameters were

p

# and the optimized ones are

p_optimized = λ2p(results.minimizer)

# After optimizing these parameters, we expect the fit to be better.
# Let's have a look at the optimized steady-state solution at the same depth as earlier, and compare it to the observed DIP field.
# First, let's calculate the steady-state with the optimized parameters.

prob_optimized = SteadyStateProblem(F, ∇ₓF, s, p_optimized)
s_optimized = solve(prob_optimized, CTKAlg()).u

# Like earlier, we fill in a 3D array of `NaN`s with the values we want to plot.
# Here, we fill the arrays with the fractional difference bewteen modeled and observed DIP:

DIPold, _ = state_to_tracers(s, nb, 2)
DIPnew, _ = state_to_tracers(s_optimized, nb, 2)
δDIPold_3D = fill(NaN, size(grd))
δDIPold_3D[iwet] .= 100(DIPold - μDIPobs) ./ μDIPobs
δDIPnew_3D = fill(NaN, size(grd))
δDIPnew_3D[iwet] .= 100(DIPnew - μDIPobs) ./ μDIPobs

# And take a slice at depth index `iz`

δDIPold_2D = δDIPold_3D[:,:,iz]
δDIPnew_2D = δDIPnew_3D[:,:,iz]

# We make those 2D slices cyclic along the longitude for Cartopy

δDIPold_cyc = hcat(δDIPold_2D, δDIPold_2D[:,1])
δDIPnew_cyc = hcat(δDIPnew_2D, δDIPnew_2D[:,1])

# And we plot the mismatch of the old, OCIM1-optimized field

figure()
δDIPlevels = -20:2:20
ax = subplot(projection = ccrs.EqualEarth(central_longitude=-155.0))
ax.coastlines()
plt2 = contourf(lon_cyc, lat, δDIPold_cyc, cmap="PiYG_r", levels=δDIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="both")
cbar2 = colorbar(plt2, orientation="horizontal", extend="both")
cbar2.set_label("δDIP / DIP [%]")
title("old DIP mismatch at $(string(round(typeof(1u"m"),grd.depth[iz]))) depth using the OCIM0.1 circulation")
gcf()

# We can see a significant negative bias in the North Atlantic.
# This is probably due to having swapped the circulation from OCIM1 to OCIM0.1.
# Let's plot the OCIM-0.1-optimized mismatch

figure()
ax = subplot(projection = ccrs.EqualEarth(central_longitude=-155.0))
ax.coastlines()
plt3 = contourf(lon_cyc, lat, δDIPnew_cyc, cmap="PiYG_r", levels=δDIPlevels, transform=ccrs.PlateCarree(), zorder=-1, extend="both")
cbar3 = colorbar(plt3, orientation="horizontal", extend="both")
cbar3.set_label("δDIP / DIP [%]")
title("new DIP mismatch at $(string(round(typeof(1u"m"),grd.depth[iz]))) depth using the OCIM0.1 circulation")
gcf()

# This shows that most of the negative bias has been corrected.
# (The green patch in the North Atlantic has dissipated.)
# Of course, this optimization was "easy", because there were only 2 tracers, few parameters, and we started from a pretty good solution.
# However, this type of optimization can really improve the skill of models, even current "state-of-the-art" ones.
# Such optimizations turn average models into better tools for cutting-edge oceanography research.

#md # !!! note
#md #     If you find any bugs or ideas you want to test and struggle with the AIBECS, don't hesitate to raise issues on the GitHub [AIBECS](https://github.com/briochemc/AIBECS.jl) repository or [contact me](www.bpasquier.com)!

#nb # > **Note**
#nb # > If you find any bugs or ideas you want to test and struggle with the AIBECS, don't hesitate to raise issues on the GitHub [AIBECS](https://github.com/briochemc/AIBECS.jl) repository or [contact me](www.bpasquier.com)!
