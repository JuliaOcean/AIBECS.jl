
# For standard use with OCIM
const iwet = (LinearIndices(wet3d))[findall(!iszero, wet3d)] # replaces find(wet3d) :(
const nwet = length(iwet)
const DIV = buildDIV(wet3d, iwet, grd)
const Iabove = buildIabove(wet3d, iwet)
const v3d = buildv3d(grd)
const v = v3d[iwet]
const vtot = sum(v)
const V = d₀(v)
vnorm²(x) = x' * V * x # Volume-weighted norm squared
Dvnorm²(x) = (2v .* x)' # Volume-weighted norm squared
vnorm(x) = sqrt(vnorm²(x))      # volume-weighted norm
vmean(x) = v'x / vtot   # volume-weighted mean
const z = grd["ZT3d"][iwet]
const ztop = grd["ZW3d"][iwet]
const spd = 24 * 60 * 60.0 # Think about using Unitful.jl

# Make up observations
const DSiobs = [10.0, 2.0, 200.0, 190.0] ;
const DSimean = vmean(DSiobs)

# Required for bgc functions
const maskEup = z .< 150.0 # Euphotic zone definition (Different from OCIM1.1!)

