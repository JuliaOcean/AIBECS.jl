

using AIBECS, Test
using SparseArrays, LinearAlgebra
using Unitful
using Unitful: °, m, km
using WorldOceanAtlasTools
using SciMLBase
using ForwardDiff, DualNumbers
using DataFrames
using Distributions
using DataDeps
using Plots

# For CI, make sure the downloads do not hang
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

test_setup_only = [:OCIM1, :OCIM0, :OCCA]
# Using `include` evaluates at global scope,
# so `Circulation` must be changed at the global scope too.
# This is why there is an `eval` in the for loop(s) below
@testset "test setup.jl only" for C in test_setup_only
    @testset "$C" begin
        eval(:(Circulation = $C))
        include("setup.jl")
    end
end


test_plots = [:OCIM2]
@testset "Test setup, plots, and sources" for C in test_plots
    @testset "$C" begin
        eval(:(Circulation = $C))
        include("setup.jl")
        include("plots.jl")
        include("sources.jl")
    end
end


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

include("AO.jl")




