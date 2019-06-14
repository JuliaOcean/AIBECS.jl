module AIBECS

using Reexport
@reexport using SparseArrays
using SuiteSparse, LinearAlgebra, Printf
using DualNumbers, DualMatrixTools
using HyperDualNumbers, HyperDualMatrixTools
using Flatten, FieldMetadata, DataFrames
using BSON
@reexport using OceanGrids            # To store the grid
@reexport using Unitful, UnitfulAstro

include("GridTools.jl")
include("SinkingParticles.jl")
include("newTypes.jl")
include("CTKsolvers.jl")
include("CTKsolvers2.jl")
include("overload_solve.jl")
include("multiTracer.jl")
include("Parameters.jl")
include("OCIM.jl")
include("SixBoxModel.jl")
include("T_2x2x2.jl")

end # module
