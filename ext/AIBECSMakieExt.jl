module AIBECSMakieExt

import AIBECS, Makie
import AIBECS: plothorizontalslice
import Unitful: ustrip, upreferred

"""

```
using AIBECS, GLMakie, JLD2
results=AIBECS.demo.demo1()
@unpack age_in_yrs,grd = results
AIBECS.plothorizontalslice(age_in_yrs,grd,depth=2000)
```
"""
function AIBECS.plothorizontalslice(x, grd; depth = nothing, kwargs...)
    tmp=AIBECS.horizontalslice(x, grd; depth = 2000);
    Makie.heatmap(ustrip(tmp))
end

end

