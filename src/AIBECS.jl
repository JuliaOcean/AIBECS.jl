module AIBECS

using Reexport
@reexport using SparseArrays
using SuiteSparse, LinearAlgebra, Printf
using ForwardDiff
@reexport using UnPack
using FieldMetadata, Flatten, DataFrames, Distributions
using BSON, NCDatasets
@reexport using OceanGrids
@reexport using Unitful
using Unitful: °
using RecipesBase, Interpolations, Distances, UnitfulRecipes
@reexport using MetadataArrays



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
include("OCIM2.jl")
include("OCCA.jl")
include("Archer_etal_2000.jl")
include("TwoBoxModel.jl")
include("Primeau_2x2x2.jl")

# AWESOME OCIM data
include("AO.jl") # TODO talk about it to Seth

# Aeolian source data
include("aeolian_sources.jl")

# Recipes for plotting
include("plot_recipes.jl")

end # module
