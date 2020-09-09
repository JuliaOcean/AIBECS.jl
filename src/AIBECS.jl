module AIBECS

using Reexport
@reexport using SparseArrays
using SuiteSparse
using LinearAlgebra
using Dates
using Printf
using ForwardDiff
@reexport using UnPack
using FieldMetadata
using Flatten
using DataFrames
using Distributions
using BSON
using NCDatasets
@reexport using OceanGrids
@reexport using Unitful
using Unitful: °
using RecipesBase
using UnitfulRecipes
using Interpolations
using Distances
using NearestNeighbors
@reexport using MetadataArrays
using Shapefile



include("CirculationGeneration.jl")
include("SinkingParticles.jl")
include("newTypes.jl")
include("CTKsolvers.jl")
include("time_steppers.jl")
include("Parameters.jl")
include("overload_solve.jl")
include("multiTracer.jl")

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

# River source data
include("Rivers.jl")

# Groundwater source data
include("GroundWaters.jl")

# ETOPO data
include("ETOPO.jl")

# Recipes for plotting
include("plot_recipes.jl")

# Diagnostics
include("diagnostics.jl")

end # module
