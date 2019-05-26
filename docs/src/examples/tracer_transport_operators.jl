
#---------------------------------------------------------
# # Tracer Transport Operators
#---------------------------------------------------------

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`tracer_transport_operators.ipynb`](@__NBVIEWER_ROOT_URL__examples/generated/tracer_transport_operators.ipynb)

# To model marine biogeochemical processes on a global scale we need to be able to account for the movement of chemical constituents both horizontally and vertically.
# We do this with a **tracer transport operator**.
# When this operator acts on a tracer field it produces the advective-diffusive divergence of the tracer.

#---------------------------------------------------------
# ## Discretization
#---------------------------------------------------------

# In order to represent the transport operator on a computer we have to discretize the tracer concentration field and the operator.
# Once discretized the tracer field is represented as a vector and the operator is represented as a sparse matrix.

#md # !!! note
#md #     A sparse matrix behaves the same way as a regular matrix.
#md #     The only difference is that in a sparse matrix the majority of the entries are zeros.
#md #     These zeros are not stored explicitly to save computer memory making it possible to deal with fairly high resolution ocean models.
#nb # > **Note**
#nb # > A sparse matrix behaves the same way as a regular matrix.
#nb # > The only difference is that in a sparse matrix the majority of the entries are zeros.
#nb # > These zeros are not stored explicitly to save computer memory making it possible to deal with fairly high resolution ocean models.

# Mathematically, the discretization converts an expression with partial derivatives into a matrix vector product:
#
# $$\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right] C \longrightarrow \mathbf{T} \, \boldsymbol{C}$$
#
# where $\mathbf{T}$ is the flux divergence transport matrix and $\boldsymbol{C}$ is the tracer concentration vector.
#
# One can go a long way towards understanding what a tracer transport operator is by playing with a simple box model. We therefore introuce a simple box model before moving on to the *Ocean Circulation Inverse Model* (OCIM).

# The simple box model we consider is embeded in a 2×2×2 "shoebox".
# It has 5 *wet* boxes and 3 *dry* boxes, as illustrated below:

#md # ```@raw html
#md # <img src="https://user-images.githubusercontent.com/4486578/58314610-3b130b80-7e53-11e9-9fe8-9527cdcca2d0.png" width =800>
#md # ```
#nb # <img src="https://user-images.githubusercontent.com/4486578/58314610-3b130b80-7e53-11e9-9fe8-9527cdcca2d0.png" width =800>

# The circulation consists of
# - a meridional overturning circulation flowing in a cycle through boxes 1 → 2 → 6 → 5 → 1 (shown in the "meridional section 1" panel above)
# - a zonal current in a reentrant cycling through boxes 1 → 3 → 1 (shown in the "layer 1" panel above)
# - vertical mixing representing deep convection between boxes 2 ↔ 6 (not shown)



#---------------------------------------------------------
# ## The transport matrix
#---------------------------------------------------------

# We start by defining the model boxes, theirs volumes, their indices, and so on:

a = 6367e3   # Earth radius           (m)
A = 4*pi*a^2 # Earth surface area     (m²)
d = 3700     # ocean depth            (m)
V = 0.75*A*d # volume of ocean        (m³)
h = 200      # thickness of top layer (m)

dz = [h*ones(4,1);(d-h)*ones(4,1)] # grid box thicknesses       (m)
dV = (dz/d).*((V/4)*ones(8,1))     # grid box volumes           (m³)
dAz = dV./dz                       # area of face ⟂ to z axis   (m²)
dy = sqrt.(dAz)                    # north-south side length    (m)
dx = sqrt.(dAz)                    # east-west side length      (m)
dAx = dV./dy                       # area of face ⟂ to x axis   (m²)
dAy = dV./dx                       # area of face ⟂ to y axis   (m²)

msk = [1, 1, 1, 0, 1, 1, 0, 0]     # wet-dry mask wet=1 dry = 0
iwet = findall(x -> x == 1, msk)        # index to wet gridboxes
idry = findall(x -> x == 0, msk)        # index to dry gridboxes
srf = [1, 1, 1, 0, 0]              # surface mask srface=1 bottom = 0
isrf = findall(x -> x == 1, srf) ;
iwet

# As you can see, `iwet` is the the vector of indices of the wet boxes.

#md # !!! note
#md #     Julia comes with Unitful, a package for using units, which AIBECS uses.
#md #     In the examples where we use AIBECS, we will use Unitful.
#nb # > **Note**
#nb # > Julia comes with Unitful, a package for using units, which AIBECS uses.
#nb # > In the examples where we use AIBECS, we will use Unitful.


# We now create the transport matrix as the flux divergence of the dissolved tracer transport due to the ocean circulation:

using LinearAlgebra
using SparseArrays
TRdiv = spzeros(8,8)
ACC = 100e6  # (m³/s) ACC for Antarctic Circumpoloar Current
TRdiv += sparse([1,1], [1,3], dV[1] \ [ACC, -ACC], 8, 8)
TRdiv += sparse([3,3], [3,1], dV[3] \ [ACC, -ACC], 8, 8)
MOC = 15e6    # (m³/s) MOC for Meridional Overturning Circulation
TRdiv += sparse([1,1], [1,5], dV[1] \ [MOC, -MOC], 8, 8)
TRdiv += sparse([2,2], [2,1], dV[2] \ [MOC, -MOC], 8, 8)
TRdiv += sparse([5,5], [5,6], dV[5] \ [MOC, -MOC], 8, 8)
TRdiv += sparse([6,6], [6,2], dV[6] \ [MOC, -MOC], 8, 8)
MIX = 10e6      # (m³/s) MIX for vertical mixing at high northern latitudes
TRdiv += sparse([2,2], [2,6], dV[2] \ [MIX, -MIX], 8, 8)
TRdiv += sparse([6,6], [6,2], dV[6] \ [MIX, -MIX], 8, 8)
TRdiv = TRdiv[iwet,iwet]
Matrix(TRdiv)

#---------------------------------------------------------
# ## Radiocarbon
#---------------------------------------------------------

# Radiocarbon, ¹⁴C, is produced by cosmic rays in the lower stratosphere and upper troposphere.
# It quickly reacts with oxygen to produce ¹⁴CO₂, which is then mixed throughout the troposphere and enters the ocean through air–sea gas exchange.
# Because the halflife of radiocarbon is only 5730 years a significant amount of decay can occur before the dissolved inorganic radiocarbon (DI¹⁴C) can mix uniformally throughout the ocean.
# As such the ¹⁴C serves as a tracer label for water that was recently in contact with the atmosphere.

# ### Tracer Equation

# Mathematically, the ¹⁴C tracer concentration, denoted $R$ (for Radiocarbon), satisfies the following tracer equation:
#
# $$\frac{\partial R}{\partial t} + \nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right] R = \Lambda(R_\mathsf{atm} - R) - R / \tau,$$
#
# where $\Lambda(R_\mathsf{atm} - R)$ represents the air–sea exchanges and $R / \tau$ the radioactive decay rate.
# ($\tau$ is the radioactive decay timescale.)
# The discretized tracer is thus given by
#
# $$\frac{\partial \boldsymbol{R}}{\partial t} + \mathbf{T} \boldsymbol{R} = \mathbf{\Lambda}(R_\mathsf{atm} - \boldsymbol{R}) - \boldsymbol{R} / \tau$$
#
# where $\mathbf{\Lambda}$ is a diagonal matrix with diagonal elements equal to zero except for surface-layer boxes.
# We can rearrange this equation as
#
# $$\frac{\partial \boldsymbol{R}}{\partial t} + \underbrace{\left[ \mathbf{T} + \mathbf{\Lambda} + \mathbf{I} / \tau \right]}_{\mathbf{M}} \, \boldsymbol{R} = \boldsymbol{s},$$
#
# where $\boldsymbol{s} = \mathbf{\Lambda} \, \boldsymbol{1} \, R_\mathsf{atm}$ effectively acts as a fixed source of Radiocarbon from the atmosphere.
# (Here we multiply $\boldsymbol{1}$, a vector of 1's, to the scalar value of $R_\mathsf{atm}$ to make sure the matrix $\mathbf{\Lambda}$ is applied to a vector, and not a scalar.)
#
# The matrix $\mathbf{M} = \mathbf{T} + \mathbf{\Lambda} + \mathbf{I} / \tau$ simplifies notation even more:
#
# $$\frac{\partial \boldsymbol{R}}{\partial t} + \mathbf{M} \, \boldsymbol{R} = \boldsymbol{s}.$$

# ### Translation to Julia Code

# We will perform an idealized radiocarbon simulation in our model and use `TRdiv`, defined earlier, for the transport matrix $\mathbf{T}$.
# In this model we prescribe the atmospheric concentration, $R_\mathsf{atm}$, to be simply equal to 1.
# (We do not specify its unit or its specific value because it is not important for determining the age of a water parcel — only the decay rate does.)
# For the air–sea gas exchange, we use a constant piston velocity $\lambda$ of 50m / 10years.

sec_per_year = 365*24*60*60 # s / yr
λ = 50 / 10sec_per_year     # m / s
Λ = λ / h * Diagonal(srf)

# where the diagonal matrix `Λ` is the discrete version of the air-sea exchange operator, $\mathbf{\Lambda}$.
# We could make it a sparse matrix via

Λ = λ / h * sparse(Diagonal(srf))

# instead, to save some memory allocations.

# For the radioactive decay we use a timescale $\tau$ of 5730/log(2) years, which we define via

τ = 5730sec_per_year / log(2)  # s

# We can now create the matrix $\mathbf{M}$ via

M = TRdiv + Λ + I / τ

# and the source term $\boldsymbol{s}$ via

s = Λ * ones(5) # air-sea source rate

# (Remember we assumed $R_\mathsf{atm}$ to be equal to 1.)


#---------------------------------------------------------
# ## Time stepping
#---------------------------------------------------------

# One way to see how the tracer evolves with time is to time step it.
# Here we will use a simple Euler-backward scheme.
# That is, we want to simulate $\boldsymbol{R}(t)$ from time $t = 0$ to time $t = \Delta t$, subject to $\boldsymbol{R}(t = 0) = 1$, using the Euler-backward scheme with $n$ timesteps.
#
# We thus discretize the $(0, \Delta t)$ time span into $n$ time intervals of size $\delta t = \Delta t / n$.
# We apply the Euler-backward scheme, which gives
#
# $$\frac{\boldsymbol{R}(t_{i+1}) - \boldsymbol{R}(t_i)}{\delta t} = - \mathbf{M} \, \boldsymbol{R}(t_{i+1}) + \boldsymbol{s}$$
#
# Reorganizing the equation to express $\boldsymbol{R}(t_{i+1})$ as a function of the rest, we get
#
# $$\big[\mathbf{I} + \delta t \, \mathbf{M}\big] \, \boldsymbol{R}(t_{i+1}) = \boldsymbol{R}(t_i) + \delta t \, \boldsymbol{s}.$$
#
# Solving such a linear system is done via the **backslash** operator, `\`.
# Let us define the time-stepping function via

euler_backward_step(R, δt, M, s) = (I + δt * M) \ (R + δt * s)

# We write another function to run all the time steps and save 

function euler_backward(R₀, ΔT, n, M, s)
    R_hist = [R₀]
    δt = Δt / n
    for i in 1:n
        push!(R_hist, euler_backward_step(R_hist[end], δt, M, s))
    end
    return reduce(hcat, R_hist), 0:δt:Δt
end

#md # !!! note
#md #     We store all the history of `R` in `R_hist`.
#md #     At each step, the `push!` function adds a new `R` to `R_hist`.
#md #     Technically, `R_hist` is a vector of vectors, so at the end, we horizontally concatenate it via `reduce(hcat, R_hist)` to rearrange it into a 2D array.
#nb # > **Note**
#nb # > We store all the history of `R` in `R_hist`.
#nb # > At each step, the `push!` function adds a new `R` to `R_hist`.
#nb # > Technically, `R_hist` is a vector of vectors, so at the end, we horizontally concatenate it via `reduce(hcat, R_hist)` to rearrange it into a 2D array.

# Now let's simulate the evolution of radiocarbon for 7500 years, starting from a concentration of 1 everywhere, via

Δt = 7500 * sec_per_year # 7500 years
R₀ = ones(5)             #
R_hist, t_hist = euler_backward(R₀, Δt, 10000, M, s) # runs the simulation

# This should take a few seconds to run.
# Once it's done, we can plot the evloution of radiocarbon through time via

using PyPlot
clf()
C14age_hist = -log.(R_hist) * τ / sec_per_year
plot(t_hist, C14age_hist')
xlabel("simulation time (years)")
ylabel("¹⁴C age (years)")
legend("box " .* string.(iwet))
title("Simulation of the evolution of ¹⁴C age with Euler-backward time steps")
gcf()

# The box model took more than 4000 years to spin up to equilibrium.
# For a box model that's no big deal because it is not computationally expensive to run the model, but for a big circulation model waiting for the model to spinup is painful.
# We therefore want a better way to find the equilibrium solution.


#---------------------------------------------------------
# ## Solving Directly for the Steady State
#---------------------------------------------------------

# One thing we notice is that when the model is at equilibrium, the $\partial \boldsymbol{R} / \partial t$ term vanishes and
# the steady state solution is simply given by the solution to the following linear system of equations
# 
# $$\mathbf{M} \, \boldsymbol{R} = \boldsymbol{s},$$
# 
# which can be solved by directly inverting the $\mathbf{M}$ matrix.

R_final = M \ s

# Converting radiocarbon into years gives the following values

C14age_final = -log.(R_final) * τ / sec_per_year
println.("box ", iwet, ": ", round.(C14age_final), " years");

# These are exactly the limit that the Radiocarbon age reaches after about 4000 years of simulation!

#---------------------------------------------------------
# ## Exercise
#---------------------------------------------------------

# Try modifying the strength of the currents of the high latitude convective mixing to see how it affects the ¹⁴C-ages.


