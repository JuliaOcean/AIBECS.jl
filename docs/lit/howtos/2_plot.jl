
#---------------------------------------------------------
# # [Plot things](@id plots)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/2_plot.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/2_plot.ipynb)

# We will create a `dummy` tracer as a function of location to showcase each plot, just for the sake of the examples herein.

# This guide is organized as follows
# - [horizontal plts](@ref horizontal-plots)
# - [vertical plots](@ref vertical-plots)
# - [Non-geographic plots](@ref non-geographic-plots)
# - 

# In this guide we will focus on how-to plot things using AIBECS' built-in recipes for [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
# These recipes are implemented using [RecipesBase.jl]https://github.com/JuliaPlots/RecipesBase.jl(), which are explained in [Plots.jl's documentation](https://docs.juliaplots.org/latest/recipes/).

# Throughout we will use the OCIM1 grid.

using AIBECS, Plots
grd, _ = OCIM1.load()

#--------------------------------------------
# ## [Horizontal plots](@id horizontal-plots)
#--------------------------------------------

# ### Horizontal slice

# The most common thing you plot after a simulation of marine tracers is a horizontal slice.
# In this case, you just need to provide the tracer (`dummy` here), the grid object `grd`, and the depth at which you want to plot.

iwet = findall(vec(iswet(grd)))
dummy = cos.(ustrip.(upreferred.(grd.lat_3D[iwet])))
horizontalslice(dummy, grd, 10)

# You can supply units for the depth at which you want to see the horizontal slice.

horizontalslice(dummy, grd, 10u"m")

# And the units should be understood under the hood.

horizontalslice(dummy, grd, 3u"km")

# If your tracer is supplied with units, those will show in the colorbar label

horizontalslice(dummy * u"mol/m^3", grd, 10u"m")

# The advantage of Plots.jl recipes like this one is that you can specify other pieces of the plot as you would with built-in functions.
# The advantage of Plots.jl recipes like this one is that you can specify other pieces of the plot as you would with built-in functions.
# For example, you can chose the colormap

dummy .= cos.(ustrip.(upreferred.(grd.lon_3D[iwet]))) .* dummy
horizontalslice(dummy, grd, 100, color=:magma)

# Or you can change predefine bits after the fact.

plt = horizontalslice(dummy, grd, 100, color=:magma)
plot!(plt, xlabel="Lon", ylabel="Lat", colorbar_title="No units", title="The pacific as a whole")

#----------------------------------------
# ## [Vertical plots](@id vertical-plots)
#----------------------------------------

# Exploring the vertical distribution of tracers is important after all.

# ### Zonal slices

# You must specify the longitude

dummy .+= sqrt.(ustrip.(grd.depth_3D[iwet])) / 30
zonalslice(dummy, grd, 330)

# ### Zonal averages

# #### Global zonal average

zonalaverage(dummy, grd)

# #### Basin zonal average

# (not available yet)

# ### Meridional slices

# (not available yet)

# ### Meridional averages

# (not available yet... at all?)

# ### Zonal transect

# (Not available yet... Thanks to GEOTRACES DMC)

# ### Meridional transect

# (Not available yet... Thanks to GEOTRACES DMC)

# ### Depth profiles

interpolateddepthprofile(dummy, grd, 0, 330)

#----------------------------------------------------
# ## [Non geographic plots](@id non-geographic-plots)
#----------------------------------------------------
