"""
    buildIbelow(wet3d, iwet)

Build the shifted-diagonal sparse matrix of the indices of below neighbours.

Ibelow[i,j] = 1 if the box represented by the linear index i
lies directly below the box represented by the linear index j.
"""
function buildIbelow(wet3d, iwet)
    nlat, nlon, ndepth = size(wet3d)
    n = nlon * nlat * (ndepth + 1)
    In = sparse(I, n, n)
    idx = zeros(Int64, nlat, nlon, ndepth + 1)
    idx[:] = 1:n
    idx .= idx[:, :, [2:ndepth + 1; 1]]      # downward shift
    return In[idx[:], :][iwet, iwet]
end

"""
    buildIabove(wet3d, iwet)

Build the shifted-diagonal sparse matrix of the indices of above neighbours.

Iabove[i,j] = 1 if the box represented by the linear index i
lies directly above the box represented by the linear index j.
"""
buildIabove(wet3d, iwet) = copy(transpose(buildIbelow(wet3d, iwet)))

buildv3d(grd) = grd["DXT3d"] .* grd["DYT3d"] .* grd["DZT3d"]

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

"""
    buildDIV(wet3d, iwet, grd)

Build the DIV operator such that

    DIV * Φ = 1/dz * (Φ - Φ[below]) ≃ dΦ/dz.
"""
function buildDIV(wet3d, iwet, grd)
    Ibelow = buildIbelow(wet3d, iwet)
    dz = grd["DZT3d"][iwet]
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
buildFLUX(w, Iabove) = sparse(Diagonal(w)) * Iabove
