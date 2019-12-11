# Reexport SinkingParticles as they are useful outside too
@reexport module SinkingParticles

using Unitful
using LinearAlgebra, SparseArrays
using ..GridTools


"""
    buildPFD(grid; settling_velocity)

Builds the particle flux divergence operator `PFD` for a given particle sinking speed
(`settling_velocity`).

Schematic of a grid cell:
```
     top  ┌─────────────────────────────────┐   ┬
          │      ↓ w                ↓ Φ     │   │
          │ (settling velovity)    (flux)   │   │
          │                                 │   │
          │            x                    │   dz
          │      (particle conc.)           │   │
          │                                 │   │
          │                                 │   │
  bottom  └─────────────────────────────────┘   ┴
```

- dz is the height of grid cell [m]
- w is the particle sinking speed at the top of the grid cell [m s⁻¹]
- Φ is the flux at the top of the grid cell [mol m⁻² s⁻¹]
- x is the particle concentration of the cell [mol m⁻³]

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
julia> buildPFD(grid, settling_velocity=1.0) # 1.0 m/s (SI units assumed)
```
"""
function buildPFD(grid; settling_velocity)
    iwet = findall(vec(grid.wet3D))
    DIV = buildDIV(grid)
    Iabove = buildIabove(grid.wet3D, iwet)
    w = ustrip.(upreferred.(settling_velocity))
    return DIV * buildFLUX(w, Iabove)
end
"""
    buildPFD(w, DIV, Iabove)

Returns `DIV * buildFLUX(w, Iabove)` 

This function is useful to avoid reconstructing `DIV` and `Iabove` every time.
It should allow for faster runs.
"""
buildPFD(w, DIV, Iabove) = DIV * buildFLUX(w, Iabove)

"""
    buildDIV(grid)

Build the DIV operator such that

    DIV * Φ = 1/dz * (Φ - Φ[below]) ≃ dΦ/dz.
"""
function buildDIV(grid)
    Ibelow = buildIbelow(grid)
    iwet = indices_of_wet_boxes(grid)
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

"""
    transportoperator(grd; kwargs)

Returns the transportoperator corresponding to the arguments.

# Example

Create the particle flux divergence with settling velocity of 100m/s
```julia-repl
julia> T = transportoperator(grd; w=100.0)
```
"""
transportoperator(grd; w, DIV=buildDIV(grd), Iabove=buildIabove(grd.wet3D, findall(vec(iswet(grd))))) = buildPFD(w, DIV, Iabove)

export transportoperator
