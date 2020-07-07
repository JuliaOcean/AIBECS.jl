```@raw html
<img src="https://user-images.githubusercontent.com/4486578/60554111-8fc27400-9d79-11e9-9ca7-6d78ee89ea70.png" alt="logo" title="AIBECS_logo" align="middle" width="50%"/>
```

# AIBECS.jl

The **Algebraic Implicit Biogeochemistry Elemental Cycling System**

!!! note
    If you are using the AIBECS for the first time, you must add it to your Julia environment, by typing `]add AIBECS` at the REPL.

This documentation is organized in 4 parts:

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

- [Parameters guide](@ref parameters)
- [Plotting basic things](@ref plots)
- [Plotting cruise/transects data](@ref cruiseplots)
- How to simulate, i.e., solve or timestep your model (coming soon!)
- How to optimize parameters (coming soon!)
- How to simulate sinking particles (coming soon!)

#### 3. Explanation/discussion

Here you will find more general discussions and explanations surrounding the AIBECS.

- [The concept of the AIBECS](@ref concept)
- [Tracer transport operators](@ref tracer-transport-operators)

#### 4. Reference

This section contains the docstrings of (almost all) the functions available in AIBECS.

----

!!! note
    The AIBECS is being developed primarily by Benoît Pasquier (postdoc with Seth John at the University of Southern California).
    If you use the AIBECS in your research, please cite it!
    [![DOI](http://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.2864051-blue.svg?&style=flat-square)](https://doi.org/10.5281/zenodo.2864051)
    Similarly, if you access data from within AIBECS (like the OCIM or OCCA ocean circulations) please cite them too.

!!! warning
    This package is in active development, so nothing is set in stone, and things may be broken sometimes.
    And if you have any suggestions or feature requests, do not hesitate to start an issue directly on the [AIBECS GitHub repository](https://github.com/JuliaOcean/AIBECS.jl), or better even, submit a pull request!
