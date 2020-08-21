```@raw html
<img src="https://user-images.githubusercontent.com/4486578/60554111-8fc27400-9d79-11e9-9ca7-6d78ee89ea70.png" alt="logo" title="AIBECS_logo" align="middle" width="50%"/>
```

# AIBECS.jl

The **Algebraic Implicit Biogeochemistry Elemental Cycling System**

!!! note
    If you are using the AIBECS for the first time, you must add it to your Julia environment, by typing
    ```
    ]add AIBECS
    ```
    at the REPL.

This documentation is organized in 4 parts:

#### 1. Tutorials

If you want to try AIBECS for the first time, this is where you should start.

- The [ideal age tutorial](@ref ideal-age) is a good place to start.
    Generate a simple linear model of an idealized tracer in a few lines of code.
- The [radiocarbon tutorial](@ref radiocarbon) is a little bit more involved, with some nonlinearities and more advanced use of AIBECS features and syntax.
- The [coupled PO₄–POP model tutorial](@ref P-model) couples 2 interacting tracers, dissolved phosphate and particulate organic phosphorus (POP).
- The [dust model tutorial](@ref dust-model) simulates some fictitious metals being injected at the ocean–atmosphere interface and being reversibly scavenged.
- The [river discharge tutorial](@ref river-discharge) similarly simulates another fictitious radioactive tracer that is injected in the ocean by rivers and that decays away as time passes.
- The [groundwater discharge tutorial](@ref groundwater-discharge) is almost identical to the river-discharge tutorial, except it uses groundwater discharge data.

#### 2. How-to guides

Here you will find goal-oriented walk-through's.

- [Parameters guide](@ref parameters)
- [Plotting basic things](@ref plots)
- [Plotting cruise/transects data](@ref cruiseplots)
- [Estimate fluxes](@ref fluxes)
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
    If you use the AIBECS in your research, [please cite it](https://doi.org/10.5281/zenodo.2864051)!
    Similarly, if you access data from within AIBECS (like the OCIM or OCCA ocean circulations) please cite them too.

!!! warning
    This package is in active development, so nothing is set in stone, and things may be broken sometimes.
    And if you have any suggestions or feature requests, do not hesitate to [start an issue](https://github.com/JuliaOcean/AIBECS.jl/issues), or better even, [submit a pull request](https://github.com/JuliaOcean/AIBECS.jl/pulls)!
