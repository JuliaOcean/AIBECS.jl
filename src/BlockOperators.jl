# Trying to come up with a lazy block type for DiffEqOpeartors
# taking mostly inspiration from the AffineDiffEqOperator in
# SciMLBase.jl/src/operators/diffeq_operator.jl



"""
HCatDiffEqOperator{T} <: AbstractBlcokDiffEqOperator{T}

```
[L1(u,p,t), L2(u,p,t)] = [Matrix(L1) Matrix(L2)] * u # after updating L1 and L2
```

"""
struct BlockDiagonalDiffEqArrayOperator{T,T1,U} <: AbstractDiffEqLinearOperator{T}
    As::T1
    du_cache::U
    function HCatDiffEqOperator{T}(As, du_cache=nothing) where T
        all([size(a) == size(As[1]) for a in As]) || error("Operator sizes do not agree")
        new{T,typeof(As),typeof(du_cache)}(As, du_cache)
    end
end

hcat(Ls::AbstractDiffEqLinearOperator...) = HCatDiffEqOperator(Ls, du_cache=Ls[1].du_cache) # not sure

Base.size(L::HCatDiffEqOperator) = (size(L.As[1],1), sum(size.(L.As, 2)))

function (L::HCatDiffEqOperator)(u,p,t::Number)
    update_coefficients!(L,u,p,t)
    du = A*u for A in L.As
    for B in L.Bs
        if typeof(B) <: Union{Number,AbstractArray}
            du .+= B
        else
            du .+= B(t)
        end
    end
    du
end

function (L::AffineDiffEqOperator)(du,u,p,t::Number)
    update_coefficients!(L,u,p,t)
    L.du_cache === nothing && error("Can only use inplace AffineDiffEqOperator if du_cache is given.")
    du_cache = L.du_cache
    fill!(du,zero(first(du)))
    # TODO: Make type-stable via recursion
    for A in L.As
        mul!(du_cache,A,u)
        du .+= du_cache
    end
    for B in L.Bs
        if typeof(B) <: Union{Number,AbstractArray}
            du .+= B
        else
            B(du_cache,t)
            du .+= du_cache
        end
    end
end

function update_coefficients!(L::AffineDiffEqOperator,u,p,t)
    # TODO: Make type-stable via recursion
    for A in L.As; update_coefficients!(A,u,p,t); end
    for B in L.Bs; update_coefficients!(B,u,p,t); end
end

@deprecate is_constant(L::AbstractDiffEqOperator) isconstant(L)
