@reexport module CirculationGeneration

using LinearAlgebra, SparseArrays
using Unitful

#=============================================
Tools for making simple circulations
=============================================#

"""
    T_advection(ϕ, orig, dest, v3D, nb)

Returns the sparse matrix transport operator, `T`.

`T` is such that it gives the flux divergence due to
the volumetric flow rate `ϕ` (in m³ s⁻¹) from box `orig` to box `dest`.
"""
function T_advection(ϕ::Tϕ, orig, dest, v3D::Array{Tv,3}, nb) where {Tϕ, Tv}
    unit(Tϕ) / unit(Tv) ≠ unit(u"1/s") && error("Unit problem")
    return sparse([orig, dest], [orig, orig], ustrip.([ϕ / v3D[orig], -ϕ / v3D[dest]]), nb, nb)
end
"""
    T_advection(ϕ, sequence, v3D, nb)

Returns the sparse matrix transport operator, `T`.

`T` is such that it gives the flux divergence due to
the volumetric flow rate `ϕ` (in m³ s⁻¹) going through all the boxes in `sequence`.
"""
function T_advection(ϕ, sequence, v3D, nb)
    origs = sequence[1:end-1]
    dests = sequence[2:end]
    sum(T_advection(ϕ, orig, dest, v3D, nb) for (orig,dest) in zip(origs, dests))
end
"""
    T_diffusion(ν, i, j, v3D, nb)

Returns the sparse matrix transport operator, `T`.

`T` is such that it gives the flux divergence due to
the volumetric mixing rate `ν` (in m³ s⁻¹) between boxes `i` and `j`.
"""
function T_diffusion(ν, i, j, v3D, nb)
    T_advection(ν, i, j, v3D, nb) + T_advection(ν, j, i, v3D, nb)
end



end
