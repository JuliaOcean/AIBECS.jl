

using AIBECS, Test
using SparseArrays, LinearAlgebra
using Unitful, UnitfulAstro
using WorldOceanAtlasTools
using DiffEqBase
using ForwardDiff
using DataFrames
# For CI, make sure the download does not hang
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

Circulation = Archer_etal_2000
@testset "Archer_etal_2000 3-box model" begin
    # Run tests with the 6-box toy model
    include("setup.jl")
end

Circulation = Primeau_2x2x2
@testset "Primeau's 6-box model" begin
    # Run tests with the 6-box toy model
    include("setup.jl")
    include("parameters.jl")
    include("bgc_functions.jl")
    include("gridtools.jl")
    include("cost_functions.jl")
    include("solvers.jl")
    include("time_steps.jl")
    include("derivatives.jl")
end

Circulation = OCIM1
@testset "OCIM1" begin
    include("setup.jl")
end

Circulation = OCIM0
@testset "OCIM0.1" begin
    include("setup.jl")
end



