#---------------------------------------------------------
# # [Parameter optimisation](@id parameter-optimization)
#---------------------------------------------------------

# This how-to describes the AIBECS interface for parameter optimisation:
# how to declare optimisable parameters, how to convert between a struct
# and an unconstrained vector, and how to evaluate the objective and its
# gradient. The end-to-end worked example — fitting a full biogeochemical
# model — lives in [GNOM](https://github.com/MTEL-USC/GNOM); this page is
# API-only.

# Optimisation features rely on package extensions that activate when you
# load the trigger packages.

using AIBECS
using Bijectors
using Distributions

# ## Declaring optimisable parameters

# Combine the `@initial_value`, `@units`, `@flattenable`, `@bounds`,
# `@logscaled`, and `@prior` macros on a single `AbstractParameters`
# struct. Each column corresponds to one piece of metadata; the
# leftmost column is consumed by the outermost macro.

# ```julia
# import AIBECS: @initial_value, @units, @flattenable, @bounds, @logscaled, @prior
# using Distributions
#
# @initial_value @units @flattenable @bounds @logscaled @prior struct OptParams{T} <: AbstractParameters{T}
#     k::T  | 1.0    | u"d^-1"  | true  | (0.01, 100.0)   | true  | LogNormal()
#     z₀::T | 100.0  | u"m"     | true  | (10.0, 5000.0)  | false | Uniform(10.0, 5000.0)
#     τ::T  | 30.0   | u"d"     | false | (1.0, 365.0)    | false | Uniform(1.0, 365.0)
# end
# ```

# Only fields with `flattenable = true` participate in optimisation; the
# others are held fixed at their initial value.

# ## Vector ↔ parameter conversion

# `vec(p)` flattens the optimisable fields to SI units. For optimisation
# we usually want an *unconstrained* vector, which is what `p2λ` produces
# by composing `vec` with each parameter's bijector (log when log-scaled,
# logit when bounded). `λ2p(typeof(p), λ)` is the inverse: it reconstructs
# a parameter struct from a flat unconstrained vector and applies inverse
# bijectors so the result satisfies the bounds.

# ```julia
# λ  = p2λ(p)
# p′ = λ2p(typeof(p), λ)
# ```

# This round-trip is exactly what `AIBECSFunction(Ts, Gs, nb, ::Type{P})`
# uses internally to let the SciML solvers receive a flat parameter
# vector even though your model is written against a struct.

# ## Objective and gradient

# The volume-weighted mismatch helper [`mismatch`](@ref) gives a single
# scalar score against an observed mean and variance (see `mismatch`):

# ```julia
# x       = ones(nb)                  # model state
# xobs    = 2.0 .* ones(nb)           # observation mean
# σ²xobs  = 0.25 .* ones(nb)          # observation variance
# v       = volumevec(grd)            # box volumes
# mismatch(x, xobs, σ²xobs, v)
# ```

# To couple the model and the data into a single objective `f(x, λ)`
# together with its analytical gradient `∇ₓf(x, λ)`, use `f_and_∇ₓf`.
# This is the entry point AIBECS uses for Newton-style
# optimisers; pair it with your favourite optimiser (e.g. via
# `Optim.jl`) to actually fit the parameters.

# ## Worked example

# A complete worked example — including building the AIBECS state
# function, hooking up `f_and_∇ₓf` to an optimiser, and reproducing
# published phosphorus + iron + dissolved-trace-metal results — is
# maintained in the GNOM repository:
#
# <https://github.com/MTEL-USC/GNOM>
