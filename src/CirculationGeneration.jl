@reexport module CirculationGeneration

using LinearAlgebra, SparseArrays
using Unitful, UnitfulAstro # for units

#=============================================
Tools for making simple circulations
=============================================#

function flux_divergence_operator_from_advection(c::Vector, i, v3D, nb)
    unit(c[1] / v3D[1]) ≠ unit(u"1/s") && error("Unit problem")
    j = circshift(i,-1) # i = origin -> j = destination
    return sparse(i, i, ustrip(c ./ v3D[i]), nb, nb) - sparse(j, i, ustrip(c ./ v3D[j]), nb, nb)
end

function flux_divergence_operator_from_advection(c, i, v3D, nb)
    c = fill(c, length(i))
    flux_divergence_operator_from_advection(c, i, v3D, nb)
end

end
