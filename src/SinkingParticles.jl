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
    return PFDO(w_top, DIVO, Iabove)
end
"""
    PFDO(w, DIVO, Iabove)

Returns `DIVO * FATO(w, Iabove)`

This function is useful to avoid reconstructing `DIVO` and `Iabove` every time.
It should allow for faster runs.
"""
PFDO(w_top, DIVO, Iabove) = DIVO * FATO(w_top, Iabove)

"""
    PFDO(grd, δz, w_top, w_bot, frac_seafloor, cumfrac_seafloor, fsedremin, Iabove)

Returns the particle-flux-divergence operator for a given sinking speed as a function of depth.

This is a slightly different construction where I take in top and bottom settling velocities,
and where the bottom one can be modified to further allow a fraction of particles to sink through
(buried into) the sea floor.

Below is a detailed explanation of how this function computes the particle-flux divergence.
Take these 3 boxes on top of each other:

┌──────────────────┐
│ water            │
│←----c----→       │
├───────────┲━━━━━━┥ ←- Φ_top
│           ┃⣿⣿⣿⣿⣿⣿│
│           ┃⣿⣿⣿⣿⣿⣿│
│←-d-→ ←-a-→┃⣿⣿⣿⣿⣿⣿│
├─────┲━━━━━┹──────┤ ←- Φ_bot
│     ┃←----b-----→│
│     ┃⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿│
│     ┃⣿⣿⣿ land ⣿⣿⣿│
│     ┃⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿│
└━━━━━┹────────────┘

A part of each of these boxes is water, while the rest is land. For the middle box,
we will use the fractional area `a = frac_seafloor` of seafloor,
and the cumulated fractional area `b = cumfrac_seafloor` of seafloor, to determine the
flux of particles at the top `ϕ_top` and bottom `ϕ_bot`.

From the top box, only the particles in the water can enter the middle box. So the flux
at the top of the middle box, `ϕ_top`, is proportional to `c = 1 - Iabove * b`.
At the bottom of the middle box, we have to take care of the case of particles hitting
the sediments and potentially buried there, or kept in the water.
For the part going through the area `d`, the flux is proportional to `d = 1 - b`.
For the part going through `a` (hitting the sediments), the flux is proportional to
`a * (1 - f)` where `f` is the fraction of particles forcibly kept in the water.
(So if `f = 0`, the particles appear to move out of the middle box, but are not carried
into the bottom one, and thus are removed from the system / buried.)
"""
function PFDO(grd, δz, w_top, w_bot, frac_seafloor, cumfrac_seafloor, fsedremin, Iabove)
    fw_bot = @. (1.0 - cumfrac_seafloor + (1.0 - fsedremin) * frac_seafloor) * w_bot
    fw_top = (1.0 .- Iabove * cumfrac_seafloor) .* w_top
    return sparse(Diagonal(fw_bot ./ δz)) - sparse(Diagonal(fw_top ./ δz)) * Iabove
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
by setting `fsedremin=0.0`.
"""
transportoperator(grd, w_top; DIVop=DIVO(grd), Iabove=buildIabove(grd)) = PFDO(w_top, DIVop, Iabove)
function transportoperator(grd, w::Function;
              δz = ustrip.(grd.δz_3D[iswet(grd)]),
              Iabove = buildIabove(grd),
              fsedremin = 1.0,
              z_top = topdepthvec(grd),
              z_bot = bottomdepthvec(grd),
              frac_seafloor = float.(isseafloorvec(grd)),
              cumfrac_seafloor = zcumsum(frac_seafloor, grd))
    return PFDO(grd, δz, ustrip.(upreferred.(w.(z_top))), ustrip.(upreferred.(w.(z_bot))), frac_seafloor, cumfrac_seafloor, fsedremin, Iabove)
end


export transportoperator

