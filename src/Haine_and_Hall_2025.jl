"""
    Haine_and_Hall_2025

9-box pedagogical circulation of Haine and Hall (2002), as described in
Figure 6 of Haine, Griffies, Gebbie, and Jiang (2025), *A review of Green's
function methods for tracer timescales and pathways in ocean models*. This
implementation translates the reference Julia version in
[OceanGreensFunctionMethods.jl](https://github.com/ggebbie/OceanGreensFunctionMethods.jl)
(file `src/pedagogical_tracer_box_model.jl`) into AIBECS's `T_advection` /
`T_diffusion` machinery. Use [`Haine_and_Hall_2025.load`](@ref) to obtain
`(grd, T)`.
"""
module Haine_and_Hall_2025
#=
Box layout: 1 longitude × 3 latitudes × 3 depths. Linear (column-major)
indexing of the wet3D = trues(3, 1, 3) array gives:

    ┌────────────────────────────────────────────┐
    │  1: H/Therm.  │  2: M/Therm.  │  3: L/Therm. │   (surface boxes 1, 2)
    ├───────────────┼───────────────┼──────────────┤
    │  4: H/Deep    │  5: M/Deep    │  6: L/Deep   │
    ├───────────────┼───────────────┼──────────────┤
    │  7: H/Abyss.  │  8: M/Abyss.  │  9: L/Abyss. │
    └────────────────────────────────────────────┘

    ╔════════════════╤════════════════╤════════════════╗
    ║                │                │                ║
    ║  1: H/Therm.   │  2: M/Therm.   │  3: L/Therm.   ║
    ║                │                │                ║
    ╟────────────────┼────────────────┼────────────────╢
    ║                │                │                ║
    ║  4: H/Deep     │  5: M/Deep     │  6: L/Deep     ║
    ║                │                │                ║
    ╟────────────────┼────────────────┼────────────────╢
    ║                │                │                ║
    ║  7: H/Abyss.   │  8: M/Abyss.   │  9: L/Abyss.   ║
    ║                │                │                ║
    ╚════════════════╧════════════════╧════════════════╝

- H = high latitudes
- M = mid latitudes
- L = low latitudes
- Therm. = thermocline
- Abyss. = abyssal
- Only boxes 1 and 2 are surface boxes
    (box 3 is considered to not exchange with the atmosphere).

Transport (all rates in Sv = 1e6 m³/s):
- abyssal overturning Ψ_a = 20 Sv on the closed loop
    3 → 2 → 1 → 4 → 7 → 8 → 9 → 6 → 3
- intermediate overturning Ψ_i = 10 Sv on the closed loop
    9 → 8 → 7 → 4 → 5 → 6 → 9
- vertical diffusive exchange ν = 5 Sv between thermocline–deep
    and deep–abyssal at every latitude.

Geometry: each box has a uniform volume of 300 Sv yr (≈ 9.5e15 m³),
matching the upstream pedagogical choice. We pick realistic depths
(0–3700 m, three equal layers) and equal-area latitude bands over
[-80°, 0°] in the southern hemisphere, then solve for the longitude
span Δlon ≈ 33° that yields the target volume.
=#

using LinearAlgebra, SparseArrays
using Unitful
using Reexport
@reexport using OceanGrids            # To store the grid
using ..CirculationGeneration
CG = CirculationGeneration


build_wet3D() = trues(3, 1, 3)


function build_grid()
    wet3D = build_wet3D()

    # Latitude: three equal-area bands over [-80°, 0°]
    Δsinlat = (sind(0) - sind(-80)) / 3
    y2 = asind(sind(0) - 1Δsinlat)
    y1 = asind(sind(0) - 2Δsinlat)
    elat = [-80.0, y1, y2, 0.0] * u"°"

    # Depth: three equal layers totalling 3700 m
    edepth = [0.0, 3700/3, 2*3700/3, 3700.0] * u"m"
    Δz = (3700/3) * u"m"

    # Longitude span chosen so each box has volume V₀ = 300 Sv·yr
    R = 6371.0u"km"
    V₀ = 300 * (1e6u"m^3/s") * 1u"yr"
    A_box = V₀ / Δz                              # required area per box
    A_full_band = 2π * R^2 * Δsinlat             # lat-band area at full longitude
    Δlon = 360u"°" * uconvert(NoUnits, A_box / A_full_band)
    elon = [0.0u"°", Δlon]

    return OceanGrid(elon, elat, edepth, wet3D)
end


function build_T(grid)
    v3D = array_of_volumes(grid)
    nb = length(v3D)

    # Abyssal overturning (Ψ_a = 20 Sv): 3 → 2 → 1 → 4 → 7 → 8 → 9 → 6 → 3
    Ψa = 20e6u"m^3/s"
    T  = CG.T_advection(Ψa, [3, 2, 1, 4, 7, 8, 9, 6, 3], v3D, nb)

    # Intermediate overturning (Ψ_i = 10 Sv): 9 → 8 → 7 → 4 → 5 → 6 → 9
    Ψi = 10e6u"m^3/s"
    T += CG.T_advection(Ψi, [9, 8, 7, 4, 5, 6, 9], v3D, nb)

    # Vertical diffusion (ν = 5 Sv) between thermocline–deep and deep–abyssal
    # at every latitude.
    ν = 5e6u"m^3/s"
    for (i, j) in [(1,4), (2,5), (3,6), (4,7), (5,8), (6,9)]
        T += CG.T_diffusion(ν, i, j, v3D, nb)
    end

    return T
end


"""
    load

Returns `grd` and `T` (in that order).
"""
function load()
    print("Creating Haine and Hall (2002) 9-box model")
    grid = build_grid()
    T = build_T(grid)
    println(" ✔")
    @info """
            You are about to use the 9-box pedagogical model of
            Haine and Hall (2002), as described in Figure 6 of the JAMES
            review paper by Haine et al. (2025).
            If you use it for research or teaching, please cite:

            - Haine, T. W. N., and T. M. Hall, 2002: A Generalized Transport Theory:
                Water-Mass Composition and Age. J. Phys. Oceanogr., 32, 1932–1946,
                https://doi.org/10.1175/1520-0485(2002)032<1932:AGTTWM>2.0.CO;2.

            - Haine, T. W. N., Griffies, S. M., Gebbie, G., & Jiang, W. (2025).
                A review of Green's function methods for tracer timescales and pathways
                in ocean models. Journal of Advances in Modeling Earth Systems, 17, e2024MS004637.
                https://doi.org/10.1029/2024MS004637

            A reference Julia implementation is at
            https://github.com/ggebbie/OceanGreensFunctionMethods.jl with the following citation:

                - Gebbie, G., Haine, T. W. N., Griffies, S. M., Jiang, W., & Pasquier, B. (2024).
                    OceanGreensFunctionMethods.jl: Julia software package with pedagogical
                    box model (v1.0). Zenodo. https://doi.org/10.5281/zenodo.14533209.
            """
    return grid, T
end

end

export Haine_and_Hall_2025
