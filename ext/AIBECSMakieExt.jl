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
plot(results)
```
"""
function plot(results::SteadyStateSolution; 
    depth = 2000, kwargs...)
    tmp=AIBECS.horizontalslice(results.state, results.grd; depth = depth);
    Makie.heatmap(ustrip(tmp))
end

end

