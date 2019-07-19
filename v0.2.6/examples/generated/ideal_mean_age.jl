using AIBECS

wet3D, grd, T_OCIM = AIBECS.OCIM1.load()
typeof(T_OCIM), size(T_OCIM)

T_age(p) = T_OCIM

source_age(age, p) = 1

z = vector_of_depths(wet3D, grd)

minimum(z)

function sink_age(age, p)
    τ = p.τ
    return age .* (z .< 20u"m") / τ
end

sms_age(age, p) = source_age(age, p) .- sink_age(age, p)

t = empty_parameter_table()    # initialize table of parameters
add_parameter!(t, :τ, 1u"s")   # add the parameter we want (τ = 1s)
initialize_Parameters_type(t, "IdealAgeParameters")  # Generate the parameter type
t

p₀ = IdealAgeParameters()

nb = number_of_wet_boxes(wet3D)  # number of wet boxes
x₀ = ones(nb)

T_matrices = (T_age,)           # bundles all the transport matrices in a tuple
sources_minus_sinks = (sms_age,) # bundles all the source-sink functions in a tuple
F, ∇ₓF = state_function_and_Jacobian(T_matrices, sources_minus_sinks, nb) # generates the state function (and its Jacobian!)
F(x₀,p₀)

prob = SteadyStateProblem(F, ∇ₓF, x₀, p₀)

age = solve(prob, CTKAlg())

iwet = indices_of_wet_boxes(wet3D)

age_3D = fill(NaN, size(wet3D)) # creates a 3D array of NaNs of the same size as `wet3D`
age_3D[iwet] = age              # Fills the wet grid boxes with the age values
size(age_3D)                    # Just to check the size of age_3D

depth = grd.depth

iz = findfirst(depth .> 1000u"m")
iz, depth[iz]

lat, lon = ustrip.(grd.lat), ustrip.(grd.lon)

age_3D_1000m_yr = age_3D[:,:,iz] * ustrip(1.0u"s" |> u"yr")

ENV["MPLBACKEND"]="qt5agg"
using PyPlot, PyCall

clf()
ccrs = pyimport("cartopy.crs")
ax = subplot(projection=ccrs.EqualEarth(central_longitude=-155.0))
ax.coastlines()
lon_cyc = [lon; 360+lon[1]] # making it cyclic for Cartopy
age_cyc = hcat(age_3D_1000m_yr, age_3D_1000m_yr[:,1])
p = contourf(lon_cyc, lat, age_cyc, levels=0:100:1200, transform=ccrs.PlateCarree(), zorder=-1)
colorbar(p, orientation="horizontal")
gcf() # gets the current figure to display

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

