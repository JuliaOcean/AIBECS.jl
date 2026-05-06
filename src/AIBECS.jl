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
@reexport using OceanGrids
@reexport using Unitful
using Unitful: °
@reexport using MetadataArrays
using SciMLBase


include("CirculationGeneration.jl")
include("SinkingParticles.jl")
include("newTypes.jl")
include("CTKsolvers.jl")
include("Parameters.jl")
include("overload_solve.jl")
include("multiTracer.jl")

# Available circulations
include("OCIM0.jl")
include("OCIM1.jl")
include("OCIM2.jl")
include("OCIM2_48L.jl")
include("OCCA.jl")
include("Archer_etal_2000.jl")
include("TwoBoxModel.jl")
include("Primeau_2x2x2.jl")
include("Haine_and_Hall_2025.jl")
export OCIM0, OCIM1, OCIM2, OCIM2_48L, OCCA
export Archer_etal_2000, TwoBoxModel
export Primeau_2x2x2, Haine_and_Hall_2025

# AWESOME OCIM data
include("AO.jl") # TODO talk about it to Seth
export AO

# Aeolian source data
include("aeolian_sources.jl")
export AeolianSources

# River source data
include("Rivers.jl")
export Rivers

# Groundwater source data
include("GroundWaters.jl")
export GroundWaters

# ETOPO data
include("ETOPO.jl")
export ETOPO

# Recipes for plotting
include("plot_recipes.jl")

# Diagnostics
include("diagnostics.jl")

end # module
