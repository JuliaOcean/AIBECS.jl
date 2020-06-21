#---------------------------------------------------------
# # [Estimate fluxes](@id fluxes)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/4_fluxes.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/4_fluxes.ipynb)


# This will take you through the process of extracting flux information from a given transport operator.
# It is split into 3 parts
# [1. Figure out the stencil of the operator](@ref operatorstencil)
# [2. Partition the operator according to said stencil](@ref operatorpartition)
# [3. Estimate the flux of a given 3D tracer field](@ref fluxestimate)

# Let's start telling Julia we will be using AIBECS and Plots.

using AIBECS, Plots

#-----------------------------------------------
# ## [Operator stencil](@id operatorstencil)
#-----------------------------------------------

# Let's load the OCIM0.1 circulation and grid

grd, T = OCIM0.load()

# and check its stencil

st = stencil(T, grd)

# This stencil shows all the indices of the neighboring cells/boxes that can exchange tracer concentrations in that model of the circulation. Let us visualize it

plotstencil(st)

# This works for any circulation, so you can swap `grd` and `T` and try again... Here we plot the same figure for a few different circulations available in AIBECS.

plts = Any[]
for Circulation in [OCIM0, OCIM1, OCIM2, OCCA]
    grd, T = Circulation.load()
    push!(plts, plotstencil(stencil(grd, T), title=string(Circulation)))
end
plot(plts..., layout=(2,2))

# These stencils are important for understanding what fluxes we are talking about.

#-----------------------------------------------
# ## [Operator Partition](@id operatorpartition)
#-----------------------------------------------

# Take a cell at index $a$ with tracer concentration $A$ and its neighbour below at index $b$ and concentration $B$.
# Depending on the concentrations $A$ and $B$, the transport operator removes (or adds) some tracer in $a$ and adds (or removes) it in $b$.
# For each neighbour (index $k$, concentration $C_k$) in the stencil, we can compute the net flow rate from $i$ to $k$ as $$M$$


# TODO
for k = 1:size(stencil,1)
    dck = stencil[k]
    if all(dck == 0) % no M_k for k=Id so skip the center stncil
        continue
    else
        iM = find(all(abs(dc - dck) == 0, 2)) ;
        M.(neighbour_name(dck)) = sparse(I(iM), J(iM), V(iM), nocn, nocn) - sparse(I(iM), I(iM), V(iM), nocn, nocn) ;
    end
end




#-----------------------------------------------
# ## [Flux Estimate](@id fluxestimate)
#-----------------------------------------------


