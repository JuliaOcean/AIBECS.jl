#---------------------------------------------------------
# # [Parameters guide](@id parameters)
#---------------------------------------------------------

# Here we will describe the AIBECS interface.
# This guide will take you through some examples of setting up model parameters.

#md # !!! note
#md #     The parameters features in AIBECS essentially come from other packages.
#md #     These include [Parameters.jl](https://github.com/mauro3/Parameters.jl), [FieldMetadata.jl](https://github.com/rafaqz/FieldMetadata.jl), [Flatten.jl](https://github.com/rafaqz/Flatten.jl), [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).

#nb # > **Note:**
#nb # > The parameters features in AIBECS essentially come from other packages.
#nb # > These include [Parameters.jl](https://github.com/mauro3/Parameters.jl), [FieldMetadata.jl](https://github.com/rafaqz/FieldMetadata.jl), [Flatten.jl](https://github.com/rafaqz/Flatten.jl), [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).

# As usual, make sure you are using AIBECS

using AIBECS
using DataFrames # required for `AIBECS.table` (loaded via package extension)

# ## Abstract parameters type

# The AIBECS provides a set of features to create and use parameters.
# This features are implemented as functions that must be provided with a parameters type.
# But each set of parameters is different, and AIBECS cannot know beforehand what parameters you want to use.
# The AIBECS thus provides an abstract parameters type, called `AbstractParameters`, upon which all the AIBECS functionality is built.
# This is why when you create a set of parameters in AIBECS, you must declare it as a subtype of `AbstractParameters`.
# Here is an example.

struct SimpleParams{T} <: AbstractParameters{T}
    α::T
    β::T
end

# Once the type, which here simply defines the symbols (`α` and `β`), is constructed, we can *instantiate* a parameter variable.

p = SimpleParams(1.0, 2.0)

# As you can see, the AIBECS will display `p` as a table.
# This is because the `show` method converts `p` to a table (a `DataFrame` to be specific) under the hood:

AIBECS.table(p)

# If you do not remember the order in which you created the parameters (`α` is first, `β` is second), AIBECS has got your back: keyword arguments are supported.

p = SimpleParams(β = 2.0, α = 1.0)

# ## Unpacking

# Probably the most useful feature in AIBECS is the ability to elegantly unpack parameters, thanks to [Parameters.jl](https://github.com/mauro3/Parameters.jl).

@unpack α, β = p
α, β

# unpacks the parameters on the left (`α` and `β`) from the parameter type `p`.

# ## Units

# One of the main features of parameters in AIBECS is that you can use units and let AIBECS do the conversions for you.
# Before you use units, though, you **must** import the `@units` and `units` functions from AIBECS.
# Here is an example.

import AIBECS: @units, units
@units struct UnitfulParams{T} <: AbstractParameters{T}
    α::T | u"m/s"
    β::T | u"d"
end

# Cretaing an instance is just as easy

p = UnitfulParams(1.0, 2.0)

# And in this case the parameters are shown with units.
# You can rely on AIBECS to convert units on the fly.

p = UnitfulParams(3.0u"km/hr", 24.0u"hr")

# And use keyword arguments

p = UnitfulParams(β = 24.0u"hr", α = 3.0u"km/hr")

# Unpacking parameters that have units first converts them to SI units.

@unpack β = p
β

# and `β` (which is equal to 1 day) here is expressed in seconds after being unpacked.

# ## Initial values

# Another useful feature is to set initial (or default) values.
# Again, you **must** import the functions for them to work properly

import AIBECS: @initial_value, initial_value
@initial_value struct ParamsWithInitialValue{T} <: AbstractParameters{T}
    α::T | 1.0
    β::T | 2.0
end

# This is handy in many applications.

# You can instantiate `p` with the initial values as, well, its values.

p = ParamsWithInitialValue()

# You could also just set one parameter to a different value

p = ParamsWithInitialValue(β = 10.0)

# ## Combining initial values and units

# You can combine both features in parameters.

@initial_value @units struct UnitfulParamsWithInitialValue{T} <: AbstractParameters{T}
    α::T | 1.0 | u"m/s"
    β::T | 2.0 | u"d"
end

# And instantiate `p` from just one parameter with its unit

p = UnitfulParamsWithInitialValue(β = 1.0u"yr")

# ## Optimizable parameters

import AIBECS: @flattenable, flattenable
@flattenable struct OptimizableParams{T} <: AbstractParameters{T}
    α::T | true
    β::T | false
    γ::T | true
end

# Then the "flattenable" parameters will be the only ones to remain when converting `OptimizableParams` to a vector

p = OptimizableParams(1.0, 2.0, 3.0)
v = vec(p)

# The `vec` function uses the `@unpack` function in AIBECS, so that units are converted when vectorizing.
# Here is an example of that by first combining `units` and `flattenable`.

@initial_value @units @flattenable struct OptimizableParamsWithUnits{T} <: AbstractParameters{T}
    α::T | 1.0 | u"m/s" | true
    β::T | 2.0 | u"d" | false
    γ::T | 3.0 | u"km" | true
end
p = OptimizableParamsWithUnits()

# And then vectorizing the parameters

vec(p)

# Note how `γ` (the third parameter, but the second flattenable one), is converted to meters.

# ## Bounds, priors, and log-scaling

# Optimisable parameters can carry per-field metadata that AIBECS uses
# during parameter optimisation: a Bayesian prior, hard bounds, and
# a flag marking the parameter as log-scaled when it spans many orders
# of magnitude. These are introduced via the `@prior`, `@bounds`, and
# `@logscaled` macros — usable on the same struct as `@units` and
# `@flattenable`. They require `using Distributions` so the
# `AIBECSParametersDistributionsExt` extension activates.
#
# ```julia
# using Distributions
# @initial_value @units @flattenable @bounds @logscaled @prior struct OptParams{T} <: AbstractParameters{T}
#     k::T | 1.0 | u"d^-1" | true  | (0.01, 100.0) | true  | LogNormal()
# end
# ```

# ## Vector ↔ parameter conversion

# Once parameters are tagged optimisable, AIBECS round-trips between a
# struct and a flat vector via `vec`, `p2λ`, and `λ2p`:
#
# - `vec(p)` flattens the optimisable fields to SI units (shown above).
# - `p2λ(p)` does the same, then maps each parameter through its bijector
#   so the vector is unconstrained (log-transformed when log-scaled,
#   logit-transformed when bounded). Pair with `λ2p(typeof(p), λ)` to
#   reconstruct a parameter struct from a flat unconstrained vector.

# ## Tabular display

# Beyond the default `show`, AIBECS exposes a few presentation helpers:
#
# - `AIBECS.table(p)` returns a `DataFrame` of the parameter values plus
#   their metadata columns (units, prior, bounds, …). Requires
#   `using DataFrames`.
# - `latex(p)` prints a LaTeX-formatted version of the same table,
#   convenient for paper supplements.
