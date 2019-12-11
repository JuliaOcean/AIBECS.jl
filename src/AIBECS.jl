module AIBECS

using Reexport
@reexport using SparseArrays
using SuiteSparse, LinearAlgebra, Printf
using ForwardDiff
@reexport using Parameters
using FieldMetadata, Flatten, DataFrames, Distributions
using BSON
@reexport using OceanGrids            # To store the grid
@reexport using Unitful, UnitfulAstro
using RecipesBase, Interpolations, Distances

include("GridTools.jl")
include("CirculationGeneration.jl")
include("SinkingParticles.jl")
include("newTypes.jl")
include("CTKsolvers.jl")
include("time_steppers.jl")
include("overload_solve.jl")
include("multiTracer.jl")
include("Parameters.jl")

# Available circulations
include("OCIM0.jl")
include("OCIM1.jl")
include("Archer_etal_2000.jl")
include("Primeau_2x2x2.jl")

# AWESOME OCIM data
include("AO.jl") # TODO talk about it to Seth

# Recipes for plotting
include("plot_recipes.jl")

end # module
