
#---------------------------------------------------------
# # [Plot basic things](@id plots)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/2_plot.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/2_plot.ipynb)


# This guide is organized as follows
# - [Horizontal maps](@ref horizontal-plots)
# - [Vertical slices](@ref vertical-plots)
# - [Depth profiles](@ref profile-plots)

# In this guide we will focus on how-to plot things using AIBECS' built-in recipes for [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
# These recipes are implemented using [RecipesBase.jl](https://github.com/JuliaPlots/RecipesBase.jl), which are explained in [Plots.jl's documentation](https://docs.juliaplots.org/latest/recipes/).

# Throughout we will use the OCIM1 grid and we will create a `dummy` tracer as a function of location to showcase each plot, just for the sake of the examples herein.

using AIBECS, Plots
grd, _ = OCIM1.load()
dummy = cosd.(latvec(grd))

#--------------------------------------------
# ## [Horizontal plots](@id horizontal-plots)
#--------------------------------------------

# ### Horizontal slice

# The most common thing you plot after a simulation of marine tracers is a horizontal slice.
# In this case, you just need to provide the tracer (`dummy` here), the grid object `grd`, and the depth at which you want to plot.

horizontalslice(dummy, grd, 10)

# You can supply units for the depth at which you want to see the horizontal slice.

horizontalslice(dummy, grd, 10u"m")

# And the units should be understood under the hood.

horizontalslice(dummy, grd, 3u"km")

# If your tracer is supplied with units, those will show in the colorbar label

horizontalslice(dummy * u"mol/m^3", grd, 10u"m")

# The advantage of Plots.jl recipes like this one is that you can specify other pieces of the plot as you would with built-in functions.
# The advantage of Plots.jl recipes like this one is that you can specify other pieces of the plot as you would with built-in functions.
# For example, you can chose the colormap with the `color` keyword argument.

dummy .*= cosd.(lonvec(grd))
plt = horizontalslice(dummy, grd, 100, color=:balance)

# And you can finetune attributes after the plot is created.

plot!(plt, xlabel="Lon", ylabel="Lat", colorbar_title="No units", title="The pacific as a whole")



#----------------------------------------
# ## [Vertical plots](@id vertical-plots)
#----------------------------------------

# Exploring the vertical distribution of tracers is important after all.

# ### Zonal slices

# You must specify the longitude

dummy .= cosd.(latvec(grd))
dummy .+= sqrt.(depthvec(grd)) / 30
zonalslice(dummy, grd, 330)

# ### Zonal averages

# #### Global zonal average

zonalaverage(dummy, grd)

# If you want a zonal average over a specific region, you can just mask it out

# #### Basin zonal average

# This is experimental at this stage and relies on an unregistered package [OceanBasins](https://github.com/briochemc/OceanBasins.jl).
# This part of the documentation will be *online* when [OceanBasins](https://github.com/briochemc/OceanBasins.jl) gets registered.

# You can create basin masks using this package with

# ```
# using OceanBasins
# mPAC = ispacific(grd)[iwet]
# surfacemap(mPAC, grd, seriestype=:heatmap, color=:lightrainbow)
# zonalaverage(dummy, grd, mask=mPAC)
# ```

# ### Meridional slices

# Just as you should expect at this stage, you can plot a meridional slice with

meridionalslice(dummy, grd, -30)

#----------------------------------------------------
# ## [Depth profiles](@id profile-plots)
#----------------------------------------------------

# Sometimes you want a profile at a given station or location

interpolateddepthprofile(dummy, grd, 0, 330)




