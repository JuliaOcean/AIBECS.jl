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

sec_per_year = 365*24*60*60
κ = 50 / 10sec_per_year     # m/s
λ  = 1 / (5730sec_per_year / log(2)); # 1/s

Λ = sparse(Diagonal(srf)) / h         # air-sea loss operator
M = TRdiv + κ * Λ + λ * I

s = κ * Λ * ones(5);             # air-sea source rate

function euler_backward(J, s, tspan, x0, n)
    dt = (tspan[2] - tspan[1]) / (n - 1)
    A = I - dt * J
    X = zeros(5, n)
    T = zeros(n)
    X[:,1] .= x0
    T[1] = tspan[1]
    for i in 2:n
        X[:,i] .= A \ (X[:,i-1] + dt * s)
        T[i] = T[i-1] + dt
    end
    return X, T
end

tspan = [0, 7500]
x0 = ones(5)
X, T = euler_backward(-sec_per_year * M, sec_per_year * s, tspan, x0, 10000)

using PyPlot
clf()
c14age = -log.(X) ./ (λ * sec_per_year)
plot(T, c14age')
xlabel("simulation time (years)")
ylabel("¹⁴C age (years)")
legend("box " .* string.(iwet))
title("Simulation of the evolution of ¹⁴C age with Euler-backward time steps")
gcf();

R = M \ s

c14age_steady_state = -log.(R) / (λ * sec_per_year)
for i in 1:5
    println("box $(iwet[i]): $(round(c14age_steady_state[i])) years")
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

