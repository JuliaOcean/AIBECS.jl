
#---------------------------------------------------------
# # [Plot transect/cruise data](@id cruiseplots)
#---------------------------------------------------------

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/3_cruiseplot.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/3_cruiseplot.ipynb)


# This guide follows the [basic plotting guide](@ref plots) and its goal is to plot data related to oceanographic observations from cruise expeditions.
# In this walkthrough, we will
# [1. Create a ficitious cruise track](@ref fictitiouscruise) with profile data
# [2. Plot transects](@ref transects)
# [3. Other plots](@ref other-plots)

# As in the [basic plotting guide](@ref plots), throughout this guide we will use the OCIM1 grid and we will create a `dummy` modelled tracer.

using AIBECS, Plots
grd, _ = OCIM1.load()
fdummy(lat, lon, depth) = @. cosd(lat) * sind(lon) + sqrt(depth) / 30
dummy = fdummy(latlondepthvecs(grd)...)

#-----------------------------------------------
# ## [A fictitious cruise](@id fictitiouscruise)
#-----------------------------------------------

# For cruise data we use the [OceanographyCruises.jl](https://github.com/briochemc/OceanographyCruises.jl) package.

#md # !!! note
#md #     The purpose of the [OceanographyCruises.jl](https://github.com/briochemc/OceanographyCruises.jl) package is to provide a Julia interface for handling discrete data from cruises.
#md #     One goal is to use it for, e.g., GEOTRACES data.
#md #     However here we merely use it to create fictitious cruise data.

#nb # > **Note:**
#nb # > The purpose of the [OceanographyCruises.jl](https://github.com/briochemc/OceanographyCruises.jl) package is to provide a Julia interface for handling discrete data from cruises.
#nb # > One goal is to use it for, e.g., GEOTRACES data.
#nb # > However here we merely use it to create fictitious cruise data.

# Let us create a station `Sydney` at (34°S, 152°E) and a `ALOHA` station at (22.75°N, 158°W).

using OceanographyCruises
Sydney = Station(name="Sydney", lat=-30, lon=156)
ALOHA = Station(name="ALOHA", lat=22.75, lon=-158)

# Now let's create a range of 10 stations from `Sydney` to `ALOHA`.

Nstations = 10
stations = range(Sydney, ALOHA, length=Nstations, westmostlon=0)

# (`westmostlon=0` ensures that the longitudes are in (0,360) to match the OCIM1 grid we use here.)

# We can now construct a fictitious cruise track

ct = CruiseTrack(name="CruisyMcCruiseFace", stations=stations)

# and check the station locations by overlaying a plot of the cruise's track over a surface map of the dummy tracer

surfacemap(dummy, grd, color=:grays)
plotcruisetrack!(ct, markercolor=:red)

# Let's create a transect of data that is almost equal to the dummy.

# First, a function for creating random depths

function randomdepths(n, max)
    depths = cumsum(rand(n+1))
    return max * view(depths,1:n) / maximum(depths)
end
Nobs = rand(1:20, Nstations) # number of obs per station/profile
depths = [randomdepths(Nobs[i], 4000) for i in 1:Nstations]
obs = [[fdummy(st.lat, st.lon, d) .+ 0.1randn() for d in depths[i]] for (i,st) in enumerate(stations)]
profiles = [DepthProfile(station=st, depths=depths[i], values=obs[i]) for (i,st) in enumerate(stations)]

t = Transect(tracer="dummy", cruise=ct.name, profiles=profiles)

#-----------------------------------------------
# ## [Transects](@id transects)
#-----------------------------------------------

# ### Zonal transect

# We can plot the modelled `dummy` data along the `ct` cruise track in the zonal directiion (along longitudes) with

zonaltransect(dummy, grd, ct)

# If we want the observations transect on top of it

zonalscattertransect!(t)

# ### Meridional transect

# Same for meridional transects (along latitude)

meridionaltransect(dummy, grd, ct)

# and

meridionalscattertransect!(t)

# If you have the GEOTRACESTools package installed and the GEOTRACES data installed at the right location, you can instead plot real data with something like
#
# ```julia
# using GEOTRACESTools
# zonalscattertransect(tracertransect("Fe", "GA02"))
# ```
#
# However, this cannot be showcased online because GEOTRACES decided its data should "not be distributed to third parties".



#----------------------------------------------------
# ## [Other plots](@id other-plots)
#----------------------------------------------------

# Plots the modelled ratio of two modelled tracers at a given station.

dummy1 = cosd.(latvec(grd)) + sqrt.(depthvec(grd)) / 30
dummy2 = cosd.(2latvec(grd)) + 0.5sqrt.(depthvec(grd)) / 30
ratioatstation(dummy1, dummy2, grd, ALOHA, xlabel="dummy 1", ylabel="dummy 2")

# This can be useful to compare stoichiometric ratios at different stations.

ratioatstation!(dummy1, dummy2, grd, Sydney, marker=:square)
