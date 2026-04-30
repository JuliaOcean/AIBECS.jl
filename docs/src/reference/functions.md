# AIBECS interface (functions and types)

## Model construction

```@docs
AIBECSFunction
split_state_function_and_Jacobian
LinearOperators
```

## Parameters and metadata

```@docs
AbstractParameters
AIBECS.table
AIBECS.latex
mismatch
f_and_∇ₓf
p2λ
λ2p
```

## Tracers and state vector

```@docs
state_to_tracers
state_to_tracer
tracers_to_state
tracer_indices
unpack_tracers
```

## Circulations

```@autodocs
Modules = [AIBECS.TwoBoxModel, AIBECS.Primeau_2x2x2, AIBECS.Archer_etal_2000, AIBECS.Haine_and_Hall_2025, AIBECS.OCIM0, AIBECS.OCIM1, AIBECS.OCIM2, AIBECS.OCIM2_48L, AIBECS.OCCA, AIBECS.CirculationGeneration]
```

## Sinking particles

```@docs
transportoperator
PFDO
DIVO
FATO
```

## Diagnostics

```@docs
stencil
directional_transport
directional_transports
smooth_operator
```

## Plot recipes

```@autodocs
Modules = [AIBECS]
Filter = t -> isa(t, Function) && (startswith(string(nameof(t)), "plot") || nameof(t) ∈ (:surfacemap, :surfacemap!, :ratioatstation, :ratioatstation!))
```

## Datasets

```@autodocs
Modules = [AIBECS.AeolianSources, AIBECS.Rivers, AIBECS.GroundWaters, AIBECS.ETOPO, AIBECS.AO]
```

## Solver

```@docs
AIBECS.CTKAlg
```

The bound `solve(::SteadyStateProblem, ::CTKAlg)` method is documented under
[`AIBECS.CTKAlg`](@ref). For the upstream `SteadyStateProblem` constructor and
generic `solve` entry point, see the
[SciMLBase documentation](https://docs.sciml.ai/SciMLBase/stable/).
