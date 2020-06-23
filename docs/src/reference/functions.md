
# AIBECS functions

## Circulations

```@docs
OCIM2.load
OCIM1.load
OCIM0.load
OCCA.load
Primeau_2x2x2.load
Archer_etal_2000.load
TwoBoxModel.load
```

## Plotting

```@docs
plothorizontalslice
surfacemap
plot∫dz
plotverticalmean
minimap
plotmeridionalslice
plotzonalmean
plot∫dx
plotzonalslice
plotmeridionalmean
plot∫dy
plot∫dxdy
plothorizontalmean
plotdepthprofile
plottransect
ratioatstation
plotparameter
plotparameters
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
```

## For optimization

```@docs
mismatch
∇mismatch
```
