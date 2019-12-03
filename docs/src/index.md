```@raw html
<img src="https://user-images.githubusercontent.com/4486578/60554111-8fc27400-9d79-11e9-9ca7-6d78ee89ea70.png" alt="logo" title="AIBECS_logo" align="middle" width="50%"/>
```

# AIBECS.jl

The **Algebraic Implicit Biogeochemistry Elemental Cycling System**

Whatever you do, if you want to use the AIBECS, you must add it to your Julia environment like every Julia package, by typing `]add AIBECS` at the REPL.

This documentation is organized in 4 parts[^1]:

[^1]: following [this documentation advice](https://www.divio.com/blog/documentation/)

#### 1. Tutorials

If you want to try AIBECS for the first time, this is where you should start.

- The [ideal age tutorial](@ref ideal-age) is a good place to start.
    It will show you how to generate a simple linear model of an idealized tracer.
- The [radiocarbon tutorial](@ref radiocarbon) is a little bit more involved,
    with some nonlinearities and more advanced use of AIBECS features and syntax
- The [coupled PO₄–POP model tutorial](@ref P-model) will show you
    how to couple 2 interacting tracers,
    one for phosphate transported by the ocean circulation,
    and one for POP transported by sinking particles.

#### 2. How-to guides

Here you will find goal-oriented walk-through's.

- How to create and use parameters in AIBECS
- [Plotting things](@ref plots)
- How to simulate, i.e., solve or timestep your model
- How to optimize parameters
- How to simulate sinking particles

#### 3. Explanation/discussion

Here you will find more general discussions and explanations surrounding the AIBECS.

- [The concept of the AIBECS](@ref concept)
- [Tracer transport operators](@ref tracer-transport-operators)

#### 4. Reference

This section contains almost all the functions available in AIBECS.


----



!!! note
    The AIBECS is being developed primarily by Benoît Pasquier with the help of François Primeau and J. Keith Moore from the Department of Earth System Science at the University of California, Irvine, and more recently with the help of Seth John from the Department of Earth Sciences at the University of Southern California.

!!! warn
    This package is in active development, so you should expect some bugs to happen. 
    And if you have any suggestions or feature requests, do not hesitate to start an issue directly on the [AIBECS GitHub repository](https://github.com/briochemc/AIBECS.jl), or even better, a submit a pull request!
