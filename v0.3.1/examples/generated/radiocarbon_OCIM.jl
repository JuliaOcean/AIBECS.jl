using AIBECS

grd, T = OCIM1.load() ;

function G(x, p)
    τ, Ratm = p.τ, p.Ratm
    return Λ(Ratm .- x, p) - x / τ
end

iwet = indices_of_wet_boxes(grd)
surface_boxes = grd.depth_3D[iwet] .== grd.depth[1]

function Λ(x, p)
    λ, h = p.λ, p.h
    return λ / h * surface_boxes .* x
end

t = empty_parameter_table()                   # initialize an empty table of parameters
add_parameter!(t, :τ, 5730u"yr"/log(2)) # radioactive decay e-folding timescale
add_parameter!(t, :λ, 50u"m" / 10u"yr") # piston velocity
add_parameter!(t, :h, grd.δdepth[1])    # height of top layer
add_parameter!(t, :Ratm, 1.0u"mol/m^3") # atmospheric concentration
t

initialize_Parameters_type(t, "C14_OCIM_parameters") # creates the type for parameters
p = C14_OCIM_parameters()                            # creates the parameters object

F, ∇ₓF = state_function_and_Jacobian(p -> T, G) # generates the state function (and its Jacobian!)
nb = length(iwet)
x = zeros(nb)
F(x,p)

prob = SteadyStateProblem(F, ∇ₓF, x, p) # define the problem
R = solve(prob, CTKAlg()).u             # solve the problem

C14age = -log.(R) * p.τ * u"s" .|> u"yr"

C14age_3D = fill(NaN, size(grd))     # creates a 3D array of NaNs

C14age_3D[iwet] .= ustrip.(C14age)   # Fills the wet grid boxes with the age values

iz = findfirst(grd.depth .> 700u"m") # aim for a depth of ~ 700 m

C14age_3D_1000m_yr = C14age_3D[:,:,iz]

ENV["MPLBACKEND"]="qt5agg"
using PyPlot, PyCall
clf()
ccrs = pyimport("cartopy.crs")
ax = subplot(projection=ccrs.EqualEarth(central_longitude=-155.0))
ax.coastlines()
lon_cyc = ustrip.([grd.lon; grd.lon[1] + 360u"°"]) # making it cyclic for Cartopy
age_cyc = hcat(C14age_3D_1000m_yr, C14age_3D_1000m_yr[:,1])
p = contourf(lon_cyc, ustrip.(grd.lat), age_cyc, levels=0:100:2000, transform=ccrs.PlateCarree(), zorder=-1)
colorbar(p, orientation="horizontal")
title("¹⁴C age at $(string(round(typeof(1u"m"),grd.depth[iz]))) depth using the OCIM1 circulation")
gcf() # gets the current figure to display

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

