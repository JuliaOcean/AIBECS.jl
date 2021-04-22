

using AIBECS, Test
using SparseArrays, LinearAlgebra
using Unitful
using Unitful: °, m, km, Myr
using WorldOceanAtlasTools
using DiffEqBase
using ForwardDiff, DualNumbers
using DataFrames
using Distributions
using Plots

# Nonlinear solver packages
using NLsolve
using NLSolvers
#using SUNDIALS # KINSOL
using SIAMFANLEquations

# For CI, make sure the downloads do not hang
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

test_setup_only = [:OCIM1, :OCIM0, :OCCA]
# Using `include` evaluates at global scope,
# so `Circulation` must be changed at the global scope too.
# This is why there is an `eval` in the for loop(s) below

println("Testing setup ---------------------------------------------------------")
@testset "test setup.jl only" for C in test_setup_only
    @testset "$C" begin
        println("  $C")
        eval(:(Circulation = $C))
        include("setup.jl")
    end
end
println("Done testing setup ----------------------------------------------------\n")


println("Testing setup, plots, and sources -------------------------------------")
test_plots = [:OCIM2]
@testset "Test setup, plots, and sources" for C in test_plots
    @testset "$C" begin
        eval(:(Circulation = $C))
        include("setup.jl")
        include("plots.jl")
        include("sources.jl")
    end
end
println("Done testing, plots, and sources --------------------------------\n")


println("Testing everything ---------------------------------------------------")
include("parameters.jl")
test_everything = [:Primeau_2x2x2, :TwoBoxModel, :Archer_etal_2000]
@testset "test everything" for C in test_everything
    @testset "$C" begin
        eval(:(Circulation = $C))
        include("setup.jl")
        include("particles.jl")
        include("bgc_functions.jl")
        include("gridtools.jl")
        include("cost_functions.jl")
        include("solvers.jl")
        include("derivatives.jl")
    end
end
println("Done testing everything ----------------------------------------------\n")






