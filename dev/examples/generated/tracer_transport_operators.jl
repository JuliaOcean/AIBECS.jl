a = 6367e3   # Earth radius           (m)
A = 4*pi*a^2 # Earth surface area     (m²)
d = 3700     # ocean depth            (m)
V = 0.75*A*d # volume of ocean        (m³)
h = 200      # thickness of top layer (m)

dz = [h*ones(4,1);(d-h)*ones(4,1)] # grid box thicknesses       (m)
dV = (dz/d).*((V/4)*ones(8,1))     # grid box volumes           (m³)
dAz = dV./dz                       # area of face ⟂ to z axis   (m²)
dy = sqrt.(dAz)                    # north-south side length    (m)
dx = sqrt.(dAz)                    # east-west side length      (m)
dAx = dV./dy                       # area of face ⟂ to x axis   (m²)
dAy = dV./dx                       # area of face ⟂ to y axis   (m²)

msk = [1, 1, 1, 0, 1, 1, 0, 0]     # wet-dry mask wet=1 dry = 0
iwet = findall(x -> x == 1, msk)        # index to wet gridboxes
idry = findall(x -> x == 0, msk)        # index to dry gridboxes
srf = [1, 1, 1, 0, 0]              # surface mask srface=1 bottom = 0
isrf = findall(x -> x == 1, srf) ;
iwet

using LinearAlgebra
using SparseArrays
TRdiv = spzeros(8,8)
ACC = 100e6  # (m³/s) ACC for Antarctic Circumpoloar Current
TRdiv += sparse([1,1], [1,3], dV[1] \ [ACC, -ACC], 8, 8)
TRdiv += sparse([3,3], [3,1], dV[3] \ [ACC, -ACC], 8, 8)
MOC = 15e6    # (m³/s) MOC for Meridional Overturning Circulation
TRdiv += sparse([1,1], [1,5], dV[1] \ [MOC, -MOC], 8, 8)
TRdiv += sparse([2,2], [2,1], dV[2] \ [MOC, -MOC], 8, 8)
TRdiv += sparse([5,5], [5,6], dV[5] \ [MOC, -MOC], 8, 8)
TRdiv += sparse([6,6], [6,2], dV[6] \ [MOC, -MOC], 8, 8)
MIX = 10e6      # (m³/s) MIX for vertical mixing at high northern latitudes
TRdiv += sparse([2,2], [2,6], dV[2] \ [MIX, -MIX], 8, 8)
TRdiv += sparse([6,6], [6,2], dV[6] \ [MIX, -MIX], 8, 8)
TRdiv = TRdiv[iwet,iwet]
Matrix(TRdiv)

sec_per_year = 365*24*60*60 # s / yr
λ = 50 / 10sec_per_year     # m / s
Λ = λ / h * Diagonal(srf)

Λ = λ / h * sparse(Diagonal(srf))

τ = 5730sec_per_year / log(2)  # s

M = TRdiv + Λ + I / τ

s = Λ * ones(5) # air-sea source rate

euler_backward_step(R, δt, M, s) = (I + δt * M) \ (R + δt * s)

function euler_backward(R₀, ΔT, n, M, s)
    R_hist = [R₀]
    δt = Δt / n
    for i in 1:n
        push!(R_hist, euler_backward_step(R_hist[end], δt, M, s))
    end
    return reduce(hcat, R_hist), 0:δt:Δt
end

Δt = 7500 * sec_per_year # 7500 years
R₀ = ones(5)             #
R_hist, t_hist = euler_backward(R₀, Δt, 10000, M, s) # runs the simulation

using PyPlot
clf()
C14age_hist = -log.(R_hist) * τ / sec_per_year
plot(t_hist, C14age_hist')
xlabel("simulation time (years)")
ylabel("¹⁴C age (years)")
legend("box " .* string.(iwet))
title("Simulation of the evolution of ¹⁴C age with Euler-backward time steps")
gcf()

R_final = M \ s

C14age_final = -log.(R_final) * τ / sec_per_year
println.("box ", iwet, ": ", round.(C14age_final), " years");

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

