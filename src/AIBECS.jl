module AIBECS

using Reexport

using WorldOceanAtlasTools
@reexport using SparseArrays
using SuiteSparse, LinearAlgebra, Printf
using DualNumbers, DualMatrixTools
using HyperDualNumbers, HyperDualMatrixTools
using Flatten, FieldMetadata, DataFrames
using JLD2
@reexport using Unitful, UnitfulAstro

include("PFDcode.jl")
include("newTypes.jl")
include("CTKsolvers.jl")
include("overload_solve.jl")
include("multiTracer.jl")
include("Parameters.jl")
include("extras.jl")
include("OCIM.jl")
include("SixBoxModelCirculation.jl")

end # module
