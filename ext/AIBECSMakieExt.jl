module AIBECSMakieExt

import AIBECS, Makie
import AIBECS.demo: SteadyStateSolution
import Unitful: ustrip, upreferred
import Makie: plot

"""
    plot(results::SteadyStateSolution; depth = 2000, kwargs...)

```
using AIBECS, GLMakie, JLD2
results=AIBECS.demo.demo1()
plot(results,depth=2000)
```
"""
function plot(results::SteadyStateSolution; 
    depth = nothing, kwargs...)
    if !isnothing(depth)
        tmp=AIBECS.horizontalslice(results.state, results.grd; depth = depth);
    else
        error("unknown option")
    end
    Makie.heatmap(ustrip(tmp))
end

end

