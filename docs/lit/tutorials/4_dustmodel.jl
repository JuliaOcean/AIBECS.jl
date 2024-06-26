
#---------------------------------------------------------
# # [A dust model](@id dust-model)
#---------------------------------------------------------

# In this tutorial we will simulate particles of dust that graviationally settle.
# 2D maps of dust deposition are available from within the AIBECS.
# To use these, we will convert them into sources in the top layer of the model grid.
# Once "born", depending on the settling velocity, these dust particles can spread horizontally until they eventually reach the sea floor and are buried in the sediments.

# The governing equation of the 3D concentration field of dust, denoted $x_\mathsf{dust}$, is:
#
# $$\left[\frac{\partial}{\partial t} + \nabla \cdot (\boldsymbol{u} + \boldsymbol{w} + \mathbf{K}\nabla)\right] x_\mathsf{dust} = s_\mathsf{dust}$$
#
# where $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \nabla \right]$ and $\nabla \cdot \boldsymbol{w}$ represent the ocean circulation transport and the vertical settling of dust particles, respectively.
# ([Tracer transport operators are described in the documentation](@ref tracer-transport-operators).)
# The only source term, $s_\mathsf{dust}$, is the dust deposition.

# In AIBECS, we must recast this equation in the generic form
#
# $$\left[\frac{\partial}{\partial t} + \mathbf{T}(\boldsymbol{p})\right] \boldsymbol{x} = \boldsymbol{G}(\boldsymbol{x},\boldsymbol{p}).$$

# We start by telling Julia we want to use the AIBECS and the OCIM0.1 transport matrix for the ocean circulation.

using AIBECS
grd, T_OCIM = OCIM0.load()

# The transport of dust is the sum of the ocean circulation and a settling-transport:

function T_dust(p)
    @unpack fsedremin = p
    return LinearOperators((T_OCIM, transportoperator(grd, z -> w(z,p), fsedremin=fsedremin)))
end

# For the settling transport of dust particles, we have used the `transportoperator` function
# for which we need to define the settling velocity `w(z,p)` as a function of depth `z` and of the parameters `p`.

#md # !!! note
#md #     We have imposed a `fsedremin` remineralization at the sediment with its keyword argument.
#md #     (Default is `fsedremin=1.0`.)
#md #     This turns the settling transport into a sink that removes all the dust reaching the sea floor.

#nb # > *Note*:
#nb # > We have imposed a `fsedremin` remineralization at the sediment with its keyword argument.
#nb # > (Default is `fsedremin=1.0`.)
#nb # > This turns the settling transport into a sink that removes all the dust reaching the sea floor.

# Following the assumption that $w(z) = w_0 + w' z$ increases linearly with depth, we write it as

function w(z,p)
    @unpack w₀, w′ = p
    return w₀ + w′ * z
end

# The parameters `fsedremin`, `w₀`, and `w′` will be defined shortly.

# ##### Dust depositon

# The AIBECS allows you to use climatological dust deposition fields from the Mahowald lab.

s_A_2D = AeolianSources.load()

# First, we must regrid these 2D maps onto the surface layer of `grd`, the 3D grid of the ocean circulation.
# To do that, we interpolate the 2D map of dust deposition onto the ocean grid's lat and lon:

s_dust_2D_monthly = s_A_2D[:dust] # kg m⁻² s⁻¹
s_dust_2D_annual = permutedims(dropdims(sum(s_dust_2D_monthly, dims=3), dims=3), (2,1)) / 12
s_dust_2D = regrid(s_dust_2D_annual, s_A_2D[:lat], s_A_2D[:lon], grd)

# (This required a bit of work to reshape `s_dust_2D_annual` before using `regrid`.)

# Then we create a blank 3D array of zeros, of which we paint the top layer:

s_dust_3D = zeros(size(grd)...)
s_dust_3D[:,:,1] .= s_dust_2D

# we convert it to a SI volumetric unit (i.e., in g m⁻³)
s_dust_3D = ustrip.(upreferred.(s_dust_3D * u"kg/m^2/s" ./ grd.δz_3D))

# Finally, we vectorize it via

s_dust = vectorize(s_dust_3D, grd) # which is the same as `s_dust_3D[iswet(grd)]`

# We then write the generic $\boldsymbol{G}$ function, which is simply the constant dust source:

G_dust(x, p) = s_dust

# ##### Parameters

# We specify some initial values for the parameters and also include units.

import AIBECS: @units, units
import AIBECS: @initial_value, initial_value
@initial_value @units struct DustModelParameters{U} <: AbstractParameters{U}
    w₀::U        | 0.1 | u"km/yr"
    w′::U        | 0.1 | u"km/yr/km"
    fsedremin::U | 5.0 | u"percent"
end

# Finally, thanks to the initial values we provided, we can instantiate the parameter vector succintly as

p = DustModelParameters()

# We build the state function `F`,

nb = count(iswet(grd))
F = AIBECSFunction(T_dust, G_dust, nb)

# the steady-state problem,

x = ones(nb) # initial guess
prob = SteadyStateProblem(F, x, p)

# and solve it

sol = solve(prob, CTKAlg()).u

# Let's now run some vizualizations using the Plots.jl recipes.

using Plots

# Let's have a look at a map of the mean dust concentration:

plotverticalmean(sol * u"g/m^3", grd, color=cgrad(:turbid, rev=true, scale=:exp))

# compared to the original source

plotverticalintegral(s_dust * u"g/m^3/s", xlabel="", ylabel="", grd, color=cgrad(:starrynight, scale=:exp), title="Surface flux", zunit=u"mg/m^2/yr")

# Huh? This is odd isn't it? What is happening? It seems that dust settles slowly enough to actually explore the ocean basins a little bit.
# This is of course dependent on our choice of model parameters.

# Let's look at what is exported below 500 m to compare:

wflux = map(z -> w(z,p), depthvec(grd)) * u"m/s"
plothorizontalslice(sol * u"g/m^3" .* wflux, grd, depth=500u"m", color=cgrad(:starrynight, scale=:exp), xlabel="", ylabel="", title="Flux at 500m", zunit=u"mg/yr/m^2")

# We can also take a look at the global mean profile (the horizontal average)

profile_plot = plothorizontalmean(sol * u"g/m^3" .|> u"mg/m^3", grd)

# As you may see, it also seems that with this model, dust (here with a 1% fraction allowed to stay in the system when it reaches the seafloor)
# accumulates a little bit in the deepest boxes of the model grid.

# Don't hesitate to play around with the parameters and see how things change!
# For example, impose a larger fraction that stays in the bottom or change other parameters via

p2 = DustModelParameters(w₀=0.1u"km/yr", w′=0.1u"km/yr/km", fsedremin=95.0u"percent")
prob2 = SteadyStateProblem(F, x, p2)
sol2 = solve(prob2, CTKAlg()).u
plotverticalmean(sol2 * u"g/m^3", grd, color=cgrad(:turbid, rev=true, scale=:exp))

# And watch how all the dust accumulates in the deepest holes of the ocean basins!
#
# Indeed, this dust tracer penetrates much deeper in the ocean:

plothorizontalmean!(profile_plot, sol2 * u"g/m^3" .|> u"mg/m^3", grd)