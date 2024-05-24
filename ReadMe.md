<a href="https://github.com/JuliaOcean/AIBECS.jl">
  <img src="https://user-images.githubusercontent.com/4486578/60554111-8fc27400-9d79-11e9-9ca7-6d78ee89ea70.png" alt="logo" title="The AIBECS logo: It represents three global marine biogeochemical cycles, where each element affects the others" align="center" width="50%"/>
</a>

# AIBECS.jl

*The ideal tool for exploring global marine biogeochemical cycles.*

<p>
  <a href="https://JuliaOcean.github.io/AIBECS.jl/stable/">
    <img src="https://img.shields.io/github/actions/workflow/status/JuliaOcean/AIBECS.jl/docs.yml?style=for-the-badge&label=Documentation&logo=Read%20the%20Docs&logoColor=white">
  </a>
  <a href="https://doi.org/10.21105/joss.03814">
    <img src="https://img.shields.io/static/v1?label=JOSS&message=10.21105/joss.03814&color=9cf&style=flat-square" alt="DOI badge">
  </a>
  <a href="https://www.bpasquier.com/talk/osm_sandiego_2020/OSM_SanDiego_2020.pdf">
    <img src=https://img.shields.io/static/v1?label=Poster&message=OSM2020&color=9cf&style=flat-square>
  </a>
</p>

<p>
  <a href="https://doi.org/10.5281/zenodo.2864051">
    <img src="http://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.2864051-blue.svg?&style=flat-square">
  </a>
  <a href="https://github.com/JuliaOcean/AIBECS.jl/blob/master/LICENSE">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-blue.svg?&style=flat-square">
  </a>
</p>

<p>
  <a href="https://github.com/JuliaOcean/AIBECS.jl/actions">
    <img src="https://img.shields.io/github/actions/workflow/status/JuliaOcean/AIBECS.jl/mac.yml?label=OSX&logo=Apple&logoColor=white&style=flat-square">
  </a>
  <a href="https://github.com/JuliaOcean/AIBECS.jl/actions">
    <img src="https://img.shields.io/github/actions/workflow/status/JuliaOcean/AIBECS.jl/linux.yml?label=Linux&logo=Linux&logoColor=white&style=flat-square">
  </a>
  <a href="https://github.com/JuliaOcean/AIBECS.jl/actions">
    <img src="https://img.shields.io/github/actions/workflow/status/JuliaOcean/AIBECS.jl/windows.yml?label=Windows&logo=Windows&logoColor=white&style=flat-square">
  </a>
  <a href="https://codecov.io/gh/JuliaOcean/AIBECS.jl">
    <img src="https://img.shields.io/codecov/c/github/JuliaOcean/AIBECS.jl/master?label=Codecov&logo=codecov&logoColor=white&style=flat-square">
  </a>
</p>




**AIBECS** (for **A**lgebraic **I**mplicit **B**iogeochemical **E**lemental **C**ycling **S**ystem, pronounced like the cool [ibex](https://en.wikipedia.org/wiki/Ibex)) is a Julia package that provides ocean biogeochemistry modellers with an easy-to-use interface for creating and running models of the ocean system.

AIBECS is a system because it allows you to choose some biogeochemical tracers, define their interactions, select an ocean circulation and *Voilà!* — your model is ready to run.

## Getting started

If you are new to AIBECS, head over to the [documentation](https://JuliaOcean.github.io/AIBECS.jl/stable/) and look for the tutorials.
(You can also click on the big "Documentation" badge above.)

## Concept

This package was developed to exploit linear-algebra tools and algorithms in Julia to efficiently simulate marine tracers.
AIBECS represents global biogeochemical cycles with a discretized system of nonlinear ordinary differential equations that takes the generic form

$$\frac{∂\boldsymbol{x}}{∂t} + \mathbf{T} \boldsymbol{x} = \boldsymbol{G}(\boldsymbol{x})$$

where $\boldsymbol{x}$ represents the model state variables, i.e., the marine tracer(s) concentration.
For a single tracer, $\boldsymbol{x}$ can be interpreted as the 3D field of its concentration.
In AIBECS, $\boldsymbol{x}$ is represented as a column vector (that's why it's **bold** and *italic*).

The operator $\mathbf{T}$ is a spatial differential operator that represents the transport of tracers.
For example, for a single tracer transported by ocean circulation,

$$\mathbf{T} = \nabla \cdot(\boldsymbol{u} - \mathbf{K}\nabla)$$

represents the effects of advection and eddy-diffusion
(where $\boldsymbol{u}$ is the 3D vector of the marine currents and $\mathbf{K}$ is a 3×3 eddy-diffusivity matrix).
Thus, $\mathbf{T}$ "acts" on $\boldsymbol{x}$ such that $\mathbf{T}\boldsymbol{x}$ is the flux divergence of that tracer.
In AIBECS, $\mathbf{T}$ is represented by a matrix (that's why it's **bold** and upstraight).

Lastly, the right-hand-side, $\boldsymbol{G}(\boldsymbol{x}$), represents the local sources minus sinks of each tracer, which must be provided as functions of the tracer(s) $\boldsymbol{x}$.

To simulate tracers using the AIBECS, you just need to define the transport operators $\mathbf{T}$ and the net sources and sinks $\boldsymbol{G}$.
That's pretty much the whole concept!

## Note on maintenance and development

Note that AIBECS is essentially in maintenance-only mode, and might remain so for a little while.
Although it works (the software has been published and has been used in publications), I consider AIBECS still in development, with a bunch of features and fixes that are currently missing, but I do not have enough time to work on it unfortunately.
Please get in touch if you have some bandwidth and/or funding and are interested in collaborating or developing AIBECS!

## References

If you use this package, please cite it.

If you use data provided by this package (like the ocean circulation from the OCIM), please cite them as well.

For convenience, all the references are available in [BibTeX](https://en.wikipedia.org/wiki/BibTeX) format in the [CITATION.bib](./CITATION.bib) file.

Also, if you want to do research using the AIBECS, and you think I could help, do not hesitate to contact me directly (contacts on my [website](www.bpasquier.com)), I would be happy to contribute and collaborate!

<img src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png" alt="NSF" title="NSF_logo" align="right" height="50"/>

The authors acknowledge funding from the Department of Energy grant DE-SC0016539 and from the National Science Foundation grant 1658380.
