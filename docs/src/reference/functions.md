
# AIBECS functions

## Circulations

```@docs
OCIM1.load
OCIM0.load
Primeau_2x2x2.load
Archer_etal_2000.load
```

## Plotting

```@docs
interpolateddepthprofile
```

## Others

```@docs
vector_of_depths
vector_of_top_depths
number_of_wet_boxes
indices_of_wet_boxes
array_of_volumes
vector_of_volumes
rearrange_into_3Darray
iswet
transportoperator
```

## Parameters

```@docs
AbstractParameters
unpack
length(<:AbstractParameters)
size(<:AbstractParameters)
values(<:AbstractParameters)
symbols
flattenable_values
flattenable_symbols
latex
table
vec
```

## For simulations

```@docs
state_function_and_Jacobian
split_state_function_and_Jacobian
SteadyStateProblem
solve
euler_forward_step
euler_forward_step!
crank_nicolson_step
crank_nicolson_step!
crank_nicolson_leapfrog_step
crank_nicolson_leapfrog_step_A⁺_and_A⁻
euler_backward_step 
euler_backward_step!
```

## For optimization

```@docs
mismatch
∇mismatch
```
