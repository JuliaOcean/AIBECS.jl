
module SixBoxModelTest

using AIBECS, Test
using SparseArrays, LinearAlgebra
using Unitful, UnitfulAstro
using WorldOceanAtlasTools
using DiffEqBase
using ForwardDiff
using DataFrames

Circulation = SixBoxModel

# For CI, make sure the download does not hang
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

@testset "Six-box model" begin
    # Run tests with the 6-box toy model
    include("setup.jl")
    include("parameters.jl")
    include("bgc_functions.jl")
    include("gridtools.jl")
    include("cost_functions.jl")
    include("solvers.jl")
    include("derivatives.jl")
end

end


