module AIBECSMakieExt

import AIBECS, Makie
import AIBECS.demo: SteadyStateSolution, @unpack
import Unitful: ustrip, upreferred
import Makie: plot

function get_axis(grd)
    lon = ustrip.(grd.lon)
    lat = ustrip.(grd.lat)
    depth = ustrip.(grd.depth)
    latticks = (-90:30:90, ["90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"])
    lonticks = (0:60:360, ["0°", "60°E", "120°E", "180°", "120°W", "60°W", "0°"])
    (lon = lon, lat = lat, depth = depth, latticks = latticks, lonticks = lonticks)
end

"""
    plot(results::SteadyStateSolution; 
        lon = nothing, lat = nothing, depth = nothing, kwargs...)

```
using AIBECS, GLMakie, JLD2
results=AIBECS.demo.demo1()
plot(results,depth=2000)
```
"""
function plot(results::SteadyStateSolution; 
    lon = nothing, lat = nothing, depth = nothing, kwargs...)

    @unpack grd, state = results
    ax=get_axis(grd)

    fig = Makie.Figure()
    if !isnothing(depth)
        slice=AIBECS.horizontalslice(state, grd; depth = depth);
        axi = Makie.Axis(fig[1, 1],title = "state at $(depth) m",
            xticks = ax.lonticks,  yticks = ax.latticks)
        hm=Makie.heatmap!(axi, ax.lon, ax.lat, slice')
    elseif !isnothing(lon)
        mer = AIBECS.meridionalslice(state, grd; lon = lon)
        axi = Makie.Axis(fig[1, 1],title = "Meridional slice at $lon",
            ylabel = "Depth (m)", xticks = ax.latticks, yreversed = true)
        hm=Makie.heatmap!(axi, ax.lat, ax.depth, mer')
    elseif !isnothing(lat)
        mer = AIBECS.zonalslice(state, grd; lat = lat)
        axi = Makie.Axis(fig[1, 1],title = "Zonal slice at $lat",
            ylabel = "Depth (m)", xticks = ax.lonticks, yreversed = true)
        hm=Makie.heatmap!(axi, ax.lon, ax.depth, mer')
    else
        error("unknown option")
    end
    Makie.Colorbar(fig[1, 2], hm)
    fig
end

end

