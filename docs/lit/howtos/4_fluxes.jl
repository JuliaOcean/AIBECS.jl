#---------------------------------------------------------
# # [Estimate fluxes](@id fluxes)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/4_fluxes.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/4_fluxes.ipynb)

# This will take you through the process of extracting flux information from a given transport operator.
# It is split into 3 parts
# 1. [Figure out the stencil of the operator](@ref operatorstencil)
# 2. [Partition the operator according to said stencil](@ref operatorpartition)
# 3. [Estimate the flux of a given 3D tracer field](@ref fluxestimate)

# Let's start telling Julia we will be using AIBECS and Plots.

using AIBECS, Plots

#-----------------------------------------------
# ## [1. Operator stencil](@id operatorstencil)
#-----------------------------------------------

# Let's load the OCCA grid information `grd` and transport matrix $\mathbf{T}$ (the variable `T` here).

grd, T = OCCA.load()

# Multiplying $\mathbf{T}$ with a given tracer concentration vector $\boldsymbol{x}$ gives the flux-divergence of that tracer.
# In other words, $\mathbf{T}$ operates on $\boldsymbol{x}$ by exchanging tracers between model-grid boxes.
# But $\mathbf{T}$ only operates on *neighboring* grid boxes.
# I.e., $\mathbf{T}$ usually does not directly exchange tracers between boxes far away from each other.

# AIBECS provides a function, `stencil`, that returns all the relative cartesian indices (AKA the relative sub-indices) of the neighbors used by a given operator.
# In the case of OCCA, the stencil is thus

st = stencil(T, grd)

# There is also a `label` function to, well, label the stencil elements in a more human-readable way.
# (Note `label` is not exported, so you must specify that it comes from the `AIBECS` package, as below.)

[st AIBECS.label.(st)]

# The stencil can be visualized with the `plotstencil` recipe

plotstencil(st)

#md # !!! note
#md #     The continuous equivalent of the matrix $\mathbf{T}$, i.e., the flux-divergence operator $\mathcal{T}$, is a differential operator, which is local by definition.
#md #     It is the spatial discretization of the differential operators that imposes the use of neighboring boxes to estimate spatial gradients.
#nb # > Note:
#nb # > In fact, the continuous equivalent of the matrix $\mathbf{T}$, i.e., the flux-divergence operator $\mathcal{T}$ is a differential operator, which is local by definition.
#nb # > It is the spatial discretization of the differential operators that imposes the use of neighboring boxes to estimate spatial gradients.

# This works for any circulation, so you can swap `grd` and `T` and try again... Here we plot the same figure for a few different circulations available in AIBECS.

plts = Any[]
for Circulation in [OCIM0, OCIM1, OCIM2, OCCA]
    grd, T = Circulation.load()
    push!(plts, plotstencil(stencil(grd, T), title=string(Circulation)))
end
plot(plts..., layout=(2,2))

# These stencils are useful to clarify the origin and destination of tracer fluxes we are going to explore in this guide..

#-----------------------------------------------
# ## [2. Operator Partition](@id operatorpartition)
#-----------------------------------------------

# For a given operator $\mathbf{T}$, to each point (or neighbor) $k$ in its stencil corresponds a unique transport operator $\mathbf{T}_k$, which effectively performs the transfer of tracers from that neighbor.
# In other words, if the stencil has $n$ elements, $\mathbf{T}$ can be partitionned into a family of $n-1$ matrices $\mathbf{T}_k$ that each operate in a single direction of the stencil.

# Let's take the "West" neighbor of the OCCA stencil

dir = st[findfirst(AIBECS.label.(st) .== "West")]

# The AIBECS provides a function `directional_T` that returns the transport operator corresponding to the direction provided from the stencil element.

T_West = directional_transport(T, grd, dir)

# Notice how `T_West` only has nonzeros on the main diagonal and on one off-diagonal, that connect each box index with the index of the box "West" of it.

#-----------------------------------------------
# ## [3. Flux Estimate](@id fluxestimate)
#-----------------------------------------------

# We can look at the depth-integrated mean flux of a fictitious tracer $\boldsymbol{x}$ of concentration 1 mol m⁻³ coming from the west by computing and plotting the vertical integral of $-\mathbf{T}\boldsymbol{x}$:

nwet = count(iswet(grd))
x = ones(nwet)
plotverticalintegral(-T_West * x * u"mol/m^3/s" .|> u"μmol/m^3/s", grd, mask=depthvec(grd) .< 100, color=:seaborn_icefire_gradient, clim=1e2 .* (-1,1))

#md # !!! warning
#md #     This guide is still a work in progress and likely contains errors
