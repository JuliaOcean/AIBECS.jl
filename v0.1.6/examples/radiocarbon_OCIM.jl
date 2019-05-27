
#---------------------------------------------------------
# # Radiocarbon with OCIM
#---------------------------------------------------------

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`radiocarbon_OCIM.ipynb`](@__NBVIEWER_ROOT_URL__examples/generated/radiocarbon_OCIM.ipynb)

using AIBECS
t = empty_parameter_table()               # initialize table of parameters
add_parameter!(t, :λ, 1 / (5730*log(2))u"yr") # add the radioactive decay e-folding timescale
add_parameter!(t, :κ, 50u"m" / 10u"yr")
initialize_Parameters_type(t, "C14_OCIM_Parameters")
t

# We generate the parameters via

p = C14_OCIM_Parameters()

# We load the OCIM matrix via

const mask, grd, T_OCIM = OCIM1.load()
T¹⁴C(p) = T_OCIM
iwet = indices_of_wet_boxes(mask);  # index to wet gridboxes

# and some useful constants

dz1 = grd["dzt"][1]               # thickness of the top layer
z = vector_of_depths(mask, grd)
nb = number_of_wet_boxes(mask)

# source minus sinks

function SMS¹⁴C(R, p) # source minus sink of age
    λ = p.λ
    κ = p.κ
    return κ * (z .< 20) .* (1.0 .- R) / dz1 - λ.*R
end

# Generate state function and Jacobian

x = ones(nb)
Ts = (T¹⁴C,)           # bundles all the transport matrices in a tuple
SMSs = (SMS¹⁴C,) # bundles all the source-sink functions in a tuple
F, ∇ₓF = state_function_and_Jacobian(Ts, SMSs, nb) # generates the state function (and its Jacobian!)
F(x,p) # to check it works

# Set up the problem in AIBECS

prob = SteadyStateProblem(F, ∇ₓF, x, p)
R = solve(prob, CTKAlg())
C14_age = -log.(R) / p.λ;

# Now we plot

lat, lon = vec(grd["yt"]), vec(grd["xt"]);
C14_age_3d = NaN * mask     # creates a 3D array of NaNs
C14_age_3d[iwet] = C14_age   # Fills the wet grid boxes with the age values
size(C14_age)            # Just to check the size

# Chose the depth index of the slice we want to plot

depth = vec(grd["zt"])
iz = findfirst(depth .> 700)
iz, depth[iz]

# get slice and convert to years

C14_age_3d_1000m_yr = C14_age_3d[:,:,iz] * ustrip(1.0u"s" |> u"yr")

# and plot

ENV["MPLBACKEND"]="qt5agg"
using PyPlot, PyCall
clf()
ccrs = pyimport("cartopy.crs")
ax = subplot(projection=ccrs.Robinson(central_longitude=-155.0))
ax.coastlines()
lon_cyc = [lon; 360+lon[1]] # making it cyclic for Cartopy
age_cyc = hcat(C14_age_3d_1000m_yr, C14_age_3d_1000m_yr[:,1])
p = contourf(lon_cyc, lat, age_cyc, levels=0:100:3600, transform=ccrs.PlateCarree(), zorder=-1)
colorbar(p, orientation="horizontal")
gcf()
