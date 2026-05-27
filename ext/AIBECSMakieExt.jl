module AIBECSMakieExt

import AIBECS, Makie
import AIBECS: plothorizontalslice
import Unitful: ustrip, upreferred

"""

```
using AIBECS, GLMakie, JLD2

grd, TOCIM2 = OCIM2.load()

function G(x, p)
    @unpack τ, z₀ = p
    return @. 1 - x / τ * (z ≤ z₀)
end

struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end

z = depthvec(grd)
nb = sum(iswet(grd))  # number of wet boxes
x_init = zeros(nb)    # Start with age = 0 everywhere

p = IdealAgeParameters(1.0, z[1])
F = AIBECSFunction(TOCIM2, G)

prob = SteadyStateProblem(F, x_init, p)
age = solve(prob, CTKAlg())
age_in_yrs = age * u"s" .|> u"yr"

AIBECS.plothorizontalslice(age_in_yrs,grd,depth=2000)
```
"""
function AIBECS.plothorizontalslice(x, grd; depth = nothing, kwargs...)
    tmp=AIBECS.horizontalslice(x, grd; depth = 2000);
    Makie.heatmap(ustrip(tmp))
end

end

