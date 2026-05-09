#---------------------------------------------------------
# # [Sinking particles](@id sinking-particles)
#---------------------------------------------------------

# This how-to walks through building transport operators for sinking
# particles — the second kind of transport in AIBECS, alongside the
# advection–diffusion of dissolved tracers. It complements the
# [Estimate fluxes](@ref fluxes) how-to, which decomposes an existing
# circulation by direction.

using AIBECS
using Plots
using JLD2 # required by `OCIM2.load`

# ## Loading a grid

# We will use OCIM2 throughout.

grd, _ = OCIM2.load()

# ## Constant settling velocity

# The simplest case: particles fall everywhere at the same speed.

T_const = transportoperator(grd, 100.0)  # 100 m/s in SI units

# `T_const` is a sparse matrix; multiplying by a tracer vector returns
# its flux divergence due to sinking.

# ## Depth-dependent settling velocity

# Pass a function `w(z)` to make the settling speed grow with depth — a
# common parameterisation. Here `w(z) = 2z + 1`, with no remineralisation
# at the seafloor (`fsedremin = 0.0` means all particles arriving at the
# seafloor are considered remineralised in place rather than reflected
# back).

T_var = transportoperator(grd, z -> 2z + 1; fsedremin = 0.0)

# ## Inspecting the building blocks

# Internally `transportoperator` is `PFDO = DIVO * FATO` — a flux-at-top
# operator (`FATO`) composed with a vertical divergence (`DIVO`):

DIVOp = DIVO(grd)
Iabove = buildIabove(grd)
F = FATO(fill(100.0, count(iswet(grd))), Iabove)
P = PFDO(fill(100.0, count(iswet(grd))), DIVOp, Iabove)

# `P` is `T_const` with the same content; exposing the pieces is useful
# when you want to compose the flux differently (e.g. attach a
# remineralisation profile, or feed `FATO` with a non-uniform `w_top`).

# ## Plotting the flux divergence

# Apply the operator to a uniform tracer field of 1 mol m⁻³ and plot a
# meridional slice of the resulting flux divergence:

x = ones(count(iswet(grd)))
divϕ = -T_var * x * u"mol/m^3/s" .|> u"μmol/m^3/d"
plotmeridionalslice(divϕ, grd, lon = 330)

# A negative divergence means the box is losing tracer to sinking out of
# the bottom; a positive divergence means it is gaining from particles
# settling in from above.
