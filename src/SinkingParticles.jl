# Reexport SinkingParticles as they are useful outside too
@reexport module SinkingParticles

using Unitful
using LinearAlgebra, SparseArrays
using ..GridTools


"""
    buildPFD(w, DIV, Iabove)

Builds the particle flux divergence operator `PFD`
for a particle sinking speed w.

Schematic of a grid cell:
```
     top  ┌─────────────────────────────┐   ┬
          │      ↓ w            ↓ Φ     │   │
          │ (sinking speed)    (flux)   │   │
          │                             │   │
          │            x                │   dz
          │      (particle conc.)       │   │
          │                             │   │
          │                             │   │
  bottom  └─────────────────────────────┘   ┴
```

- dz is the height of grid cell [m]
- w is the particle sinking speed at the top of the grid cell [m s⁻¹]
- Φ is the flux at the top of the grid cell [mmol m⁻² s⁻¹]
- x is the particle concentration of the cell [mmol m⁻³]

The PFD operator is defined by

    PFD * x ≃ dΦ/dz,

i.e., so that applied to x it approximates dΦ/dz.
It is calculated as the matrix product of the DIV and FLUX operators:

    PFD = DIV * FLUX.

DIV, the divergence operator (d/dz only here), is defined by

    DIV * Φ = 1/dz * (Φ - Φ[below]) ≃ dΦ/dz.

And FLUX, the sinking particle flux operator, is defined by

    FLUX * x = w * x[above] ≃ Φ.

# Example

```julia-repl
julia> PFD = buildPFD(w)
```
"""
buildPFD(w, DIV, Iabove) = DIV * buildFLUX(w, Iabove)

function buildPFD(grid, wet3D; sinking_speed=1.0)
    w = sinking_speed
    iwet = findall(vec(wet3D))
    DIV = buildDIV(wet3D, iwet, grid)
    Iabove = buildIabove(wet3D, iwet)
    return DIV * buildFLUX(w, Iabove)
end

"""
    buildDIV(wet3D, iwet, grid)

Build the DIV operator such that

    DIV * Φ = 1/dz * (Φ - Φ[below]) ≃ dΦ/dz.
"""
function buildDIV(wet3D, iwet, grid)
    Ibelow = buildIbelow(wet3D, iwet)
    dz = ustrip.(grid.δz_3D[iwet])
    return sparse(Diagonal(dz.^(-1))) * (Ibelow - I) # divergence with positive downwards
end

"""
    buildFLUX(w, Iabove)

Build the `FLUX` operator for a particle sinking speed ``w(r)``

Note `w` is the sinking speed at the top of each grid cell.
(`w` can also be a constant (Float64 type).)

The `FLUX` operator is defined by

    FLUX * x = w * x(above) ≃ Φ.
"""
buildFLUX(w::Vector, Iabove) = sparse(Diagonal(w)) * Iabove
buildFLUX(w::Number, Iabove) = w * Iabove

export buildDIV, buildPFD, buildFLUX

end
