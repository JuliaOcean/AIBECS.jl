# AIBECS.jl Documentation

(In construction)

This package provides XXX.


## Usage

```@meta
DocTestSetup = quote
    using AIBECS
    using LinearAlgebra, DiffEqBase
end
```

Load the OCIM 1.1 matrix and grid with 

```jldoctest usage
wet3d, grd, T_OCIM = AIBECS.OCIM1.load()
typeof(T_OCIM)

# output

SparseMatrixCSC{Float64,Int64}
```


