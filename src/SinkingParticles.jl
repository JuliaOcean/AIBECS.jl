# Reexport SinkingParticles as they are useful outside too
#@reexport module SinkingParticles

using Unitful
using LinearAlgebra, SparseArrays
using OceanGrids


"""
    PFDO(grd; w_top)

Builds the particle-flux-divergence operator `PFDO` for a given particle sinking speed
(`w_top`).

Schematic of a grid cell:
```
     top  ┌─────────────────────────────────┐   ┬
          │      ↓ w_top          ↓ Φ_top   │   │
          │ (settling velovity)    (flux)   │   │
          │                                 │   │
          │            x                    │   δz
          │      (particle conc.)           │   │
          │                                 │   │
          │                                 │   │
  bottom  └─────────────────────────────────┘   ┴
```

- `δz` is the height of grid cell [m]
- `w_top` is the particle sinking speed at the top of the grid cell [m s⁻¹]
- `Φ_top` is the flux at the top of the grid cell [mol m⁻² s⁻¹]
- `x` is the particle concentration of the cell [mol m⁻³]

The PFDO is defined by

    PFDO * x = δΦ/δz ≈ dΦ/dz,

i.e., so that applied to `x` it approximates the flux divergence of `x`, dΦ/δz.
It is calculated as the matrix product of DIVO and FATO:

    PFDO = DIVO * FATO.

where the divergence operator `DIVO`, is defined by (`z` increasing downward)

    DIVO * ϕ_top = 1/δz * (ϕ_top[below] - ϕ_top) = δϕ/δz ≈ dΦ/δz.

and FATO, the flux-at-top operator, is defined by

    FATO * x = w_top * x[above] ≈ ϕ_top.

# Example

```julia-repl
julia> PFDO(grd, w_top=1.0) # 1.0 m/s (SI units assumed)
```
"""
function PFDO(grd; w_top)
    iwet = findall(vec(grd.wet3D))
    DIVO = DIVO(grd)
    Iabove = buildIabove(grd.wet3D, iwet)
    w_top = ustrip.(upreferred.(w_top))
    return DIVO * FATO(w_top, Iabove)
end
"""
    PFDO(w, DIVO, Iabove)

Returns `DIVO * FATO(w, Iabove)`

This function is useful to avoid reconstructing `DIVO` and `Iabove` every time.
It should allow for faster runs.
"""
PFDO(w_top, DIVO, Iabove) = DIVO * FATO(w_top, Iabove)

"""
    PFDO(grd, w::Function)

Returns the particle-flux-divergence operator for a given sinking speed as a function of depth.

This is a slightly different construction where I allow `w` to be used as a function of depth,
which further allows particles to sink through (buried into) the sea floor.
"""
function PFDO(grd, w::Function;
              δz = ustrip.(grd.δz_3D[iswet(grd)]), Iabove=buildIabove(grd), sedremin=true,
              w_top = ustrip.(upreferred.(w.(topdepthvec(grd)))),
              w_bot = (sedremin ? Iabove' * ones(length(w_bot)) : 1) .* ustrip.(upreferred.(w.(bottomdepthvec(grd)))))
    return sparse(Diagonal(w_bot ./ δz)) - sparse(Diagonal(w_top ./ δz)) * Iabove
end


"""
    DIVO(grd)

Build the `DIVO` operator such that

    DIVO * ϕ_top = 1/δz * (ϕ_top - ϕ_top[below]) ≈ dΦ/δz.
"""
function DIVO(grd)
    Ibelow = buildIbelow(grd)
    iwet = indices_of_wet_boxes(grd)
    δz = ustrip.(grd.δz_3D[iwet])
    return sparse(Diagonal(1 ./ δz)) * (Ibelow - I) # divergence with positive downwards
end

"""
    FATO(w_top, Iabove)

Build the `FATO` operator for a particle sinking speed `w_top`

(`w_top` is the sinking speed at the top of each grid cell.)

The `FATO` operator is defined by

    FATO * x = w_top * x(above) ≈ ϕ_top.
"""
FATO(w_top::Vector, Iabove) = sparse(Diagonal(w_top)) * Iabove
FATO(w_top::Number, Iabove) = w_top * Iabove

export DIVO, PFDO, FATO
#end # module

"""
    transportoperator(grd; kwargs)

Returns the transportoperator corresponding to the arguments.

# Examples

Create the particle flux divergence with settling velocity of 100m/s

```julia-repl
julia> T = transportoperator(grd; w = 100.0)
```

Or with settling velocity function w(z) = 2z + 1

```julia-repl
julia> T = transportoperator(grd; w = z -> 2z + 1)
```

By default, the seafloor flux is set to zero, so that all the particles
that reach it are remineralized there. You can let particles go through
by setting `sedremin=false`.
"""
transportoperator(grd, w_top; DIVop=DIVO(grd), Iabove=buildIabove(grd)) = PFDO(w_top, DIVop, Iabove)
function transportoperator(grd, w::Function; sedremin=true, z_top=topdepthvec(grd), z_bot=bottomdepthvec(grd), 
                           Iabove=buildIabove(grd), δz=ustrip.(grd.δz_3D[iswet(grd)]))
    return PFDO(grd, w; δz=δz, Iabove=Iabove, sedremin=sedremin, 
                w_top = ustrip.(upreferred.(w.(z_top))),
                w_bot = (sedremin ? Iabove' * ones(length(z_bot)) : 1) .* ustrip.(upreferred.(w.(z_bot))))
end

export transportoperator

