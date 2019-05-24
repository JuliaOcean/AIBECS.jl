
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
# ## From Equation to Code
#---------------------------------------------------------

# Radiocarbon, ¹⁴C, is produced by cosmic rays in the lower stratosphere and upper troposphere.
# It quickly reacts with oxygen to produce ¹⁴CO₂, which is then mixed throughout the troposphere and enters the ocean through air-sea gas exchange.
# Because the halflife of radiocarbon is only 5730 years a significant amount of deday can occur before the dissolved inorganic radiocarbon (DI¹⁴C) can mix uniformally throughout the ocean.
# As such the ¹⁴C serves as a tracer label for water that was recently in contact with the atmosphere.

# ### Tracer Equation

# Mathematically, the ¹⁴C tracer concentration, denoted $R$ (for Radiocarbon), satisfies the following tracer equation:
#
# $$\frac{\partial R}{\partial t} = - \nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right] R + \Lambda(R - [^{14}C]_\mathsf{atm}) - \kappa R$$
#
# The discretized tracer is thus given by
#
# $$\frac{\partial \boldsymbol{R}}{\partial t} = - \mathbf{T} \boldsymbol{R} + \Lambda(\boldsymbol{R} - [^{14}C]_\mathsf{atm}) - \kappa \boldsymbol{R}$$
#
# We can rearrange that equation as
#
# $$\frac{dR}{dt} + \left[\mathbf{T}+\lambda\mathbf{I}+\kappa\boldsymbol{\Lambda}\right]R = \boldsymbol{s},$$
#
# where $\boldsymbol{s} = \kappa\boldsymbol{\Lambda}\boldsymbol{1}$ effectively acts as a fixed source of Radiocarbon from the atmosphere.


# ### Translation to Julia Code

# Here we will perform an idealized radiocarbon simulation in our model and use `TRdiv` for the transport matrix $\mathbf{T}$.
# In this model we prescribe the atmospheric concentration to 1 (unit?) and model the air-sea gas exchange using a constant piston velocity $\kappa$ of 50m / 10years.
# For the radioactive decay we use a timescale $\tau$ of (5730 years) / log(2).

sec_per_year = 365*24*60*60
κ = 50 / 10sec_per_year     # m/s
λ  = 1 / (5730sec_per_year / log(2)); # 1/s

#-

Λ = sparse(Diagonal(srf)) / h         # air-sea loss operator
M = TRdiv + κ * Λ + λ * I

#- 

s = κ * Λ * ones(5);             # air-sea source rate




#---------------------------------------------------------
# ## Simulating Radiocarbon
#---------------------------------------------------------

# One way to see how the tracer evolves with time is to time step it.
# Here we will use a simple Euler-backward scheme.
# That is, below we create a function to solve $dx/dt = F(x)$ from `tspan[1]` to `tspan[2]` subject to `x[tspan[1]] = x0` using the Euler-backward scheme with `n` timesteps.
# Note that `(x[i+1]-x[i])/dt = J*X[i+1] + s`.

function euler_backward(J, s, tspan, x0, n)
    dt = (tspan[2] - tspan[1]) / (n - 1)
    A = I - dt * J
    X = zeros(5, n)
    T = zeros(n)
    X[:,1] .= x0
    T[1] = tspan[1]
    for i in 2:n
        X[:,i] .= A \ (X[:,i-1] + dt * s)
        T[i] = T[i-1] + dt
    end
    return X, T
end

#-

tspan = [0, 7500]
x0 = ones(5)
X, T = euler_backward(-sec_per_year * M, sec_per_year * s, tspan, x0, 10000)

#-

using PyPlot
clf()
c14age = -log.(X) ./ (λ * sec_per_year)
plot(T, c14age')
xlabel("simulation time (years)")
ylabel("¹⁴C age (years)")
legend("box " .* string.(iwet))
title("Simulation of the evolution of ¹⁴C age with Euler-backward time steps")
gcf();

# The box model took more than 4000 years to spin up to equilibrium.
# For a box model that's no big deal because it is not computationally expensive to run the model, but for a big circulation model waiting for the model to spinup is painful.
# We therefore want a better way to find the equilibrium solution.


#---------------------------------------------------------
# ## Solving for the Steady State
#---------------------------------------------------------

# One thing we notice is that when the model is at equilibrium, the $dR/dt$ term vanishes and
# the steady state solution is given by the solution to the following linear system of equations
#
# $$\underbrace{\left[\mathbf{T}+\lambda\mathbf{I}+\kappa\boldsymbol{\Lambda}\right]}_{\mathbf{M}}R = \underbrace{\kappa\boldsymbol{\Lambda}\mathbf{1}}_{s}$$
#
# which can be solved by directly inverting the $M$ matrix.

R = M \ s

#-

c14age_steady_state = -log.(R) / (λ * sec_per_year)
for i in 1:5
    println("box $(iwet[i]): $(round(c14age_steady_state[i])) years")
end

# These are exactly the limit that the Radiocarbon age reaches after about 4000 years of simulation!
# Try modifying the strength of the currents of the high latitude convective mixing to see how it affects the ¹⁴C-ages.


