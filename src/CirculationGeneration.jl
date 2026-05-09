@reexport module CirculationGeneration

using LinearAlgebra
using SparseArrays
using Unitful

"""
    T_advection(ϕ, orig, dest, v3D, nb)

Returns the sparse matrix transport operator, `T`.

`T` is such that it gives the flux divergence due to
the volumetric flow rate `ϕ` (in m³ s⁻¹) from box `orig` to box `dest`.
"""
function T_advection(ϕ::Tϕ, orig, dest, v3D::Array{Tv, 3}, nb) where {Tϕ, Tv}
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
    origs = sequence[1:(end - 1)]
    dests = sequence[2:end]
    sum(T_advection(ϕ, orig, dest, v3D, nb) for (orig, dest) in zip(origs, dests))
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

@doc """
    CirculationGeneration

Tools for assembling simple advection/diffusion transport operators from
sequences of inter-box flow rates, with `T_advection` and `T_diffusion`
as the two building blocks. Used by the pedagogical circulation modules
([`TwoBoxModel`](@ref), [`Primeau_2x2x2`](@ref), [`Archer_etal_2000`](@ref),
[`Haine_and_Hall_2025`](@ref)).
""" CirculationGeneration
