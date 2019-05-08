using AIBECS, Test

# 4-box model tests
include("build_4box_model_circulation.jl") # define constants T, grd, and wet3d
@test size(T, 1) == 4
@test size(T, 2) == 4

include("setup_4box_model.jl") # define constants T, grd, and wet3d
@test vmean(ones(nwet)) ≈ 1.0
@test vnorm²(ones(nwet)) ≈ 1.292e18     # m³   # From Archer et al. [2000]


# # Load wet3d, grid and T (OCIM type)
# # const path_to_OCIM = "/Users/UCIben/Dropbox/OCIM1.1/CTL.mat"
# # const wet3d, grd, T = loadOCIM(path_to_OCIM)
# using JLD
# # @load "/Users/UCIben/Google Drive/JuliaData/OCIM1.1/CTL.jld" # loads wet3d grd T
# # @save "/Users/UCIben/Google Drive/JuliaData/OCIM1.1/CTL.jld" wet3d grd T
# # @save "/Users/UCIben/.julia/v0.6/OCIMtools/data/OCIM1.1_CTL.jld" wet3d grd T
# @load (homedir() * "/.julia/v0.6/OCIMtools/data/OCIM1.1_CTL.jld") wet3d grd T

# const iwet = find(wet3d)
# const nwet = length(iwet)

# # Test the constructing functions

# const DIV = buildDIV(wet3d, grd)
# const Iabove = buildIabove(wet3d)
# const spd = 24 * 60 * 60.0

# const τg = 1e6 * 365 * spd
# const v3d = buildv3d(grd)
# const v = v3d[iwet]
# const vtot = sum(v)
# const vn = v / vtot
# const svn = sqrt.(vn)

# # include(homedir() * "/.julia/v0.6/OCIMtools/test/transportMatrixTests.jl")

# # const path_to_obs = "/Users/UCIben/Projects/CSD2/results/data/DSi_obs_OCIMgrid.mat"
# # const obs_name = "DSi_obs"
# # const xobs = load_obs(path_to_obs,obs_name)
# # @save "/Users/UCIben/.julia/v0.6/OCIMtools/data/OCIM1.1_SiObs.jld" xobs
# @load (homedir() * "/.julia/v0.6/OCIMtools/data/OCIM1.1_SiObs.jld") xobs
# const xgeo = vn ⋅ xobs
# geores(x) = (xgeo - x) / τg

# # using Parameters.jl for my parameters
# using Parameters.jl

# @with_kw struct Para
#     τu = 30.0 * spd # specific uptake rate timescale, default = 30 d
#     w₀ = 10 / spd   # terminal sinking velocity at surface, default = 10 m/d
#     Dzw = 1 / spd   # dw/dz, default = 1 (m/d)/m
# end

# function paraprint2(p::Para)
#     @unpack τu, w₀, Dzw = p
#     @printf "Parameter values\n"
#     @printf "  τu = %.2g days\n" τu / spd
#     @printf "  w₀ = %.2g m d⁻¹\n" w₀ * spd
#     @printf "  ∂w/∂z = %.2g m d⁻¹ / m\n" Dzw * spd
# end
# # end of Parameters.jl trial

# τu(p) = p[1]
# w₀(p) = p[2]
# Dzw(p) = p[3]
# function paraprint(p)
#     @printf "Parameter values\n"
#     @printf "  τu = %.2g days\n" τu(p) / spd
#     @printf "  w₀ = %.2g m d⁻¹\n" w₀(p) * spd
#     @printf "  ∂w/∂z = %.2g m d⁻¹ / m\n" Dzw(p) * spd
# end
# function paraprint(p::Array{Complex{Float64}})
#     paraprint(real.(p))
#     print("Imaginary part: ")
#     paraprint(imag.(p))
# end
# using DualNumbers
# function paraprint(p::Array{DualNumbers.Dual{Float64}})
#     paraprint(realpart.(p))
#     print("Epsilon part: ")
#     paraprint(epsilon.(p))
# end

# const p₀ = [30 * spd, 10 / spd, 1 / spd]
# const x₀ = xgeo * ones(nwet)
# const np = length(p₀)

# include(homedir() * "/.julia/v0.6/OCIMtools/test/derivativeTests.jl")

# const d0 = spdiagm
# const τ = 30.0 * spd
# const z = grd["ZT3d"][iwet]
# const ztop = grd["ZW3d"][iwet]
# const ieup = z .< 85.0

# xpos(x) = (x .>= 0) .* x
# xpos(x::Array{Complex{Float64},1}) = (real.(x) .>= 0.0) .* x
# xpos(x::Array{DualNumbers.Dual{Float64},1}) = (realpart.(x) .>= 0.0) .* x
# Dxpos(x) = d0((x .>= 0) .* 1.0)
# Dxpos(x::Array{Complex{Float64},1}) = d0((real.(x) .>= 0) .* 1.0)
# Dxpos(x::Array{DualNumbers.Dual{Float64},1}) = d0((realpart.(x) .>= 0) .* 1.0)
# derivTest(x -> xpos(x), x -> Dxpos(x), 1e2 * (2rand(nwet) - 1))

# ubase(x, p) = ieup .* xpos(x - xobs) ./ τu(p)
# Dxubase(x, p) = d0(ieup) * Dxpos(x - xobs) ./ τu(p)
# derivTest(x -> ubase(x, p₀), x -> Dxubase(x, p₀), 1e2 * (2rand(nwet) - 1))

# u(x, p) = ubase(xpos(x), p)
# Dxu(x, p) = Dxubase(xpos(x), p) * Dxpos(x)
# derivTest(x -> u(x, p₀), x -> Dxu(x, p₀), 1e2 * (2rand(nwet) - 1))

# function Dpu(x, p)
#     Dpbuf = zeros(eltype(p), nwet, np)
#     Dpbuf[:,1] .= -u(x, p) ./ τu(p)
#     Dpbuf[:,2:3] .= 0
#     return Dpbuf
# end
# derivTest(p -> u(x₀, p), p -> Dpu(x₀, p), p₀)

# # const w = 100 / spd # m s⁻¹ = m d⁻¹ /(s d⁻¹)
# # const S = buildPFD(w,DIV,Iabove)


# const i_c, j_c, ij_c, i_d, j_d, ij_d = buildCompressionIndices(wet3d)
# comp(X) = compress(X, i_c, j_c, ij_d)
# decomp(Xc) = decompress(Xc, i_d, j_d, ij_c)

# const κ = 0.25 / spd
# const S1 = buildPFD(ones(nwet), DIV, Iabove)
# const Sz = buildPFD(ztop, DIV, Iabove)
# S(p) = buildPFD(w₀(p) + Dzw(p) * ztop, DIV, Iabove)

# mutable struct SFactors
#     Sf
#     p
# end
# Sfc = SFactors(factorize(S(p₀) + κ * I), p₀)
# function update_Sf(p, Sfc)
#     if (eltype(p) != eltype(Sfc.p)) || (p != Sfc.p)
#         Sfc.p = p
#         Sfc.Sf = factorize(S(p) + κ * I)
#     end
#     return Sfc
# end
# function particlex(x, p)
#     update_Sf(p, Sfc)
#     return Sfc.Sf \ u(x, p)
# end
# function Dxparticlex(x, p)
#     update_Sf(p, Sfc)
#     return decomp(Sfc.Sf \ comp(Dxu(x, p)))
# end
# function Dpparticlex(x, p)
#     update_Sf(p, Sfc)
#     Dpbuf = Sfc.Sf \ Dpu(x, p)
#     Dpbuf[:,2] .= Sfc.Sf \ (S1 * -particlex(x, p))
#     Dpbuf[:,3] .= Sfc.Sf \ (Sz * -particlex(x, p))
#     return Dpbuf
# end
# complexStepDerivTest(x -> particlex(x, p₀), x -> Dxparticlex(x, p₀), x₀)
# complexStepDerivTest(p -> particlex(x₀, p), p -> Dpparticlex(x₀, p), p₀)
# derivTest(x -> particlex(x, p₀), x -> Dxparticlex(x, p₀), x₀)
# derivTest(p -> particlex(x₀, p), p -> Dpparticlex(x₀, p), p₀)

# r(x, p) = κ * particlex(x, p)
# Dxr(x, p) = κ * Dxparticlex(x, p)
# Dpr(x, p) = κ * Dpparticlex(x, p)
# complexStepDerivTest(x -> r(x, p₀), x -> Dxr(x, p₀), x₀)
# complexStepDerivTest(p -> r(x₀, p), p -> Dpr(x₀, p), p₀)
# derivTest(x -> r(x, p₀), x -> Dxr(x, p₀), x₀)
# derivTest(p -> r(x₀, p), p -> Dpr(x₀, p), p₀)

# # state function
# f(x, p) = -T * x + r(x, p) - u(x, p) + geores(x)
# # TODO think about making use of u(x,p) only once
# Dxf(x, p) = -T + (Dxr(x, p) - Dxu(x, p)) - I / τg
# complexStepDerivTest(x -> f(x, p₀), x -> Dxf(x, p₀), 1e2 * (2rand(nwet) - 1))
# derivTest(x -> f(x, p₀), x -> Dxf(x, p₀), 1e2 * (2rand(nwet) - 1))

# Dpf(x, p) = Dpr(x, p) - Dpu(x, p)
# complexStepDerivTest(p -> f(x₀, p), p -> Dpf(x₀, p), p₀)
# derivTest(p -> f(x₀, p), p -> Dpf(x₀, p), p₀)

# # TODO write proper mass conservation check
# v ⋅ (r(x₀, p₀) - u(x₀, p₀)) / (v ⋅ x₀)
# v ⋅ (T * x₀) / (v ⋅ x₀)
# v ⋅ (r(x₀, p₀) - u(x₀, p₀)) / (v ⋅ abs.(r(x₀, p₀) - u(x₀, p₀)))

# # Sanity checks
# # FDxf(x,p) = sparse(ForwardDiff.jacobian(x -> f(x,p),x))
# # TODO write ForwardDiff Jacobian
# #FDpf(x,p) = sparse(ForwardDiff.jacobian(p -> f(x,p),p))

# # cost functions
# δ(x) = (x - xobs) / xgeo
# Dδ(x) = 1 / xgeo
# complexStepDerivTest(δ, Dδ, x₀)
# derivTest(δ, Dδ, x₀)

# # cost function
# c(x) = 0.5 * vn ⋅ δ(x).^2
# Dc(x) = (vn .* δ(x) .* Dδ(x)).'
# complexStepDerivTest(c, Dc, x₀)
# foo, bar = complexStepDerivDebug(c, Dc, x₀)
# abs((foo - bar) / bar)
# derivTest(c, Dc, x₀)
# foo, bar = derivDebug(c, Dc, x₀)
# abs((foo - bar) / bar)
# #
# # # Sanity check
# # using ForwardDiff
# # FDc(x) = ForwardDiff.gradient(c,x).'

# nrm(x) = norm(svn .* x)
# const τstop = τg * 1e6

# # Real test begins here
# x₀ .= xgeo * ones(nwet)
# SaJ = StateAndJacobian(x₀, factorize(Dxf(x₀, p₀)), p₀) # the Jacobian factors

# function qprint!(c, p, SaJ, f, Dxf, nrm, τstop)
#     paraprint(p)
#     return q!(c, p, SaJ, f, Dxf, nrm, τstop)
# end
# qwrap(p) = qprint!(c, p, SaJ, f, Dxf, nrm, τstop)
# # qwrap(p) = q!(c, p, SaJ, f, Dxf, nrm, τstop)
# # qwrap(p₀)

# xsol = solvef0NewtonChordShamanskii!(x -> f(x, p₀), x -> Dxf(x, p₀), x₀, nrm, τstop, true)

# complexStepDerivTest(p -> f(xsol, p), p -> Dpf(xsol, p), p₀)
# complexStepDerivTest(x -> f(x, p₀), x -> Dxf(x, p₀), xsol)
# complexStepDerivTest(c, Dc, xsol)
# derivTest(p -> f(xsol, p), p -> Dpf(xsol, p), p₀)
# derivTest(x -> f(x, p₀), x -> Dxf(x, p₀), xsol)
# derivTest(c, Dc, xsol)

# # x₀ .= xgeo * ones(nwet)
# # f3d = NaN * wet3d
# # f3d[iwet] .= u(x₀,p₀)
# # using Plots
# # gr()
# # x₀
# # sum(u(x₀,p₀).<0) / sum(u(x₀,p₀).>0)
# #
# # plot(f3d[30,50:60,:]')
# #
# #
# # xsol = ΨtcNewtonSolve!(x->f(x,p₀), x->Dxf(x,p₀), x₀, nrm, τstop)



# const λ₀ = log.(p₀)
# # Need to define the function with a storage argument first
# function Qprint!(c, λ, SaJ, f, Dxf, nrm, τstop)
#     paraprint(exp.(λ))
#     qval = Q!(c, λ, SaJ, f, Dxf, nrm, τstop)
#     @printf "=====> RMS = %g%%\n\n" sqrt(2qval)
#     return qval
# end
# Qwrap(λ) = Qprint!(c, λ, SaJ, f, Dxf, nrm, τstop)
# # Qwrap(λ) = Q!(c, λ, SaJ, f, Dxf, nrm, τstop)
# slowQwrap(λ) = slowQ(c, λ, nwet, f, Dxf, nrm, τstop)
# # nrm(f(ones(length(x₀)),p₀))
# # Q₀ = Qwrap(λ₀)
# DQwrap(λ) = DQ!(Dc, λ, SaJ, f, Dxf, Dpf, nrm, τstop)
# slowDQwrap(λ) = slowDQ(Dc, λ, nwet, f, Dxf, Dpf, nrm, τstop)
# D2Qwrap(λ) = D2Q!(Dc, λ, SaJ, f, Dxf, Dpf, nrm, τstop)
# # @test derivTest(Qwrap, DQwrap, λ₀)
# # @test complexStepDerivTest(slowQwrap, slowDQwrap, λ₀)

# # Δλ = sqrt.(sqrt.(eps.(λ₀))) .* bitrand(np)
# # Qwrap(λ₀ + Δλ) - Qwrap(λ₀)
# # DQwrap(λ₀) * Δλ
# #
# # print_with_color(:red,"slowq  test:\n")
# # Dqwrap(p) = Dq!(Dc, p, SaJ, f, Dxf, Dpf, nrm, τstop)
# # Δp = sqrt(eps()) * p₀ .* bitrand(np)
# # p₀
# # p₀ + Δp
# # (p₀ + Δp) - p₀
# # print_with_color(:red,"  q2:\n")
# # q2 = qwrap(p₀ + Δp)
# # print_with_color(:red,"  q1:\n")
# # q1 = qwrap(p₀)
# # q2-q1
# # (q2-q1)/q1
# # print_with_color(:red,"  Dq * Δp:\n")
# # Dq1 = Dqwrap(p₀)
# # Dq1 * Δp


# # print_with_color(:red,"New test:\n")
# # slowqwrap(p) = slowq(c, p, nwet, f, Dxf, nrm, τstop)
# # slowDqwrap(p) = slowDq(Dc, p, nwet, f, Dxf, Dpf, nrm, τstop)
# # Δp = sqrt(eps()) * p₀ .* bitrand(np)
# # p₀
# # p₀ + Δp
# # (p₀ + Δp) - p₀
# # print_with_color(:red,"  q2:\n")
# # q2 = slowqwrap(p₀ + Δp)
# # print_with_color(:red,"  q1:\n")
# # q1 = slowqwrap(p₀)
# # q2-q1
# # (q2-q1)/q1
# # print_with_color(:red,"  Dq * Δp:\n")
# # Dq1 = slowDqwrap(p₀)
# # Dq1 * Δp

# # print_with_color(:red,"CSD test:\n")
# # Δp1 = im * h * [1,0]
# # Δp2 = im * h * [0,1]
# # CSDq1 = imag.(slowqwrap(p₀ + Δp1)) / h
# # CSDq2 = imag.(slowqwrap(p₀ + Δp2)) / h
# # Dq1 = slowDqwrap(p₀)

# # print_with_color(:red,"CSD fast test:\n")
# # Δp1 = im * h * [1,0]
# # Δp2 = im * h * [0,1]
# # CSDq1 = imag.(qwrap(p₀ + Δp1)) / h
# # CSDq2 = imag.(qwrap(p₀ + Δp2)) / h
# # Dq1 = Dqwrap(p₀)




# # complexStepDerivTest(Qwrap, DQwrap, λ₀)
# # SaJ.x

# function DQstor!(storage, λ)
#     storage[1:np] .= DQwrap(λ).'
# end
# function D2Qstor!(storage, λ)
#     storage[1:np, 1:np] .= D2Qwrap(λ)
# end


# function print_results(results)
#     λopt = results.minimizer
#     popt = exp.(λopt)
#     paraprint(popt)
#     qopt = qwrap(popt)
#     print_with_color(:red, "  gives cost q = $(@sprintf("%.2g",qopt))\n")
#   #println("    which means ≈ $(@sprintf("%.2g",(100*sqrt(qopt))))\% RMS")
#   #println("      And steady state a∞ = $(newton_solver(popt))")
#   #println(results)
# end

# using Optim
# # print("\nOptimizing (Q Newton)...")
# # @benchmark print_results(optimize(Qwrap, λ₀, Newton()))
# # print("\nOptimizing (Newton with FD gradient and Hessian)...")
# # @benchmark print_results(optimize(slowQ, FDQ!, FD2Q!, λ₀, Newton()))
# # print("\nOptimizing (Newton with gradient and Hessian)...")
# # print_results(optimize(Qwrap, DQstor!, D2Qstor!, λ₀, Newton()))

# x₀ .= xgeo / 100 * ones(nwet)
# @printf "\n\nTesting newton split with debug\n"
# xsol = solvef0NewtonChordShamanskii!(x -> f(x, p₀), x -> Dxf(x, p₀), x₀, nrm, τstop, true)

# x₀ .= xgeo / 100 * ones(nwet)
# @printf "\n\nTesting pseudo-transient coninutatioin\n"
# xsol = ψtcNewtonSolve!(x -> f(x, p₀), x -> Dxf(x, p₀), x₀, nrm, τstop)




# # # Benchmarks
# # using BenchmarkTools
# #
# # @benchmark slowQwrap(λ₀)
# #
# # @benchmark print_results(optimize(Qwrap, DQstor!, D2Qstor!, λ₀, NewtonTrustRegion()))
# # #
# # @benchmark print_results(optimize(Qwrap, λ₀, BFGS()))
