using AIBECS
mask, grd, T_OCIM = OCIM1.load() ;

iwet = findall(vec(mask))       # index to wet grdboxes
h = ustrip(grd.δdepth[1])     # thickness of the top layer
z = ustrip.(grd.depth_3D[iwet]) # depth of the grdbox centers
nwet = length(iwet)             # number of wet grd boxes
v = vector_of_volumes(mask,grd)# volume of the grdboxes
lat, lon = ustrip.(grd.lat), ustrip.(grd.lon) # latitudes and longitudes (useful for plotting)
depth = ustrip.(grd.depth) ;

t = empty_parameter_table()               # initialize table of parameters
add_parameter!(t, :λ, 1 / (5730*log(2))u"yr") # add the radioactive decay e-folding timescale
add_parameter!(t, :κ, 50u"m" / 10u"yr")
add_parameter!(t, :z₀, 20u"m")
initialize_Parameters_type(t, "C14_Parameters") # Generate the parameter table

p = C14_Parameters()

function sms(R, p)
    λ, κ, z₀ = p.λ, p.κ, p.z₀
    return κ / h * (z .< z₀) .* (1.0 .- R) - λ * R
end

T_14c(p) = T_OCIM

T_matrices = (T_14c,)            # bundles all the transport matrices in a tuple
sources_minus_sinks = (sms,) # bundles all the source-sink functions in a tuple
F, ∇ₓF = state_function_and_Jacobian(T_matrices, sources_minus_sinks, nwet) # generates the state function (and its Jacobian!)

x = ones(nwet)                           # initial iterate for the solver
prob = SteadyStateProblem(F, ∇ₓF, x, p) # define the problem
R = solve(prob, CTKAlg())                 # solve the problem
c14age = -log.(R) / p.λ ;                   # convert R to $^{14}C-age$

c14age_3D = fill(NaN, size(mask))     # creates a 3D array of NaNs
c14age_3D[iwet] = c14age   # Fills the wet grd boxes with the age values
size(c14age_3D)            # Just to check the size of age_3D

iz = findfirst(depth .> 700) # aim for a depth of ~ 700 m
depth_map = Int(round(depth[iz]))
iz, depth_map
c14age_3D_1000m_yr = c14age_3D[:,:,iz] * ustrip(1.0u"s" |> u"yr")

ENV["MPLBACKEND"]="qt5agg"
using PyPlot, PyCall

clf()
ccrs = pyimport("cartopy.crs")
ax = subplot(projection=ccrs.EqualEarth(central_longitude=-155.0))
ax.coastlines()
lon_cyc = [lon; 360+lon[1]] # making it cyclic for Cartopy
age_cyc = hcat(c14age_3D_1000m_yr, c14age_3D_1000m_yr[:,1])
p = contourf(lon_cyc, lat, age_cyc, levels=0:100:1800, transform=ccrs.PlateCarree(), zorder=-1)
colorbar(p, orientation="horizontal")
title("¹⁴C age at $(depth_map)m depth")
gcf() # gets the current figure to display

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

