<a href="https://github.com/briochemc/AIBECS.jl">
  <img src="https://user-images.githubusercontent.com/4486578/60554111-8fc27400-9d79-11e9-9ca7-6d78ee89ea70.png" alt="logo" title="The AIBECS logo: It represents three global marine biogeochemical cycles, where each element affects the others" align="center" width="50%"/>
</a>

# AIBECS.jl

*The ideal tool for exploring global marine biogeochemical cycles.*

<p>
  <a href="https://briochemc.github.io/AIBECS.jl/stable/">
    <img src=https://img.shields.io/badge/docs-stable-important.svg?style=flat-square&label=Documentation&logo=Read%20the%20Docs>
  </a>
</p>

<p>
  <a href="https://doi.org/10.5281/zenodo.2864051">
    <img src="http://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.2864051-blue.svg?&style=flat-square">
  </a>
  <a href="https://github.com/briochemc/AIBECS.jl/blob/master/LICENSE">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-blue.svg?&style=flat-square">
  </a>
</p>

<p>
  <a href="https://travis-ci.com/briochemc/AIBECS.jl">
    <img alt="Build Status" src="https://img.shields.io/travis/com/briochemc/AIBECS.jl/master?label=OSX/Linux/Windows&logo=travis&logocolor=white&style=flat-square">
  </a>
  <a href="https://codecov.io/gh/briochemc/AIBECS.jl">
    <img src="https://img.shields.io/codecov/c/github/briochemc/AIBECS.jl/master?label=Codecov&logo=codecov&logoColor=white&style=flat-square">
  </a>
</p>






**AIBECS** (for **A**lgebraic **I**mplicit **B**iogeochemical **E**lemental **C**ycling **S**ystem, pronounced like the cool [ibex](https://en.wikipedia.org/wiki/Ibex)) is a Julia package that provides ocean biogeochemistry modelers with an easy-to-use interface for creating and running models of the ocean system.

AIBECS is a system because it allows you to chose some biogeochemical tracers, define their interactions, select an ocean circulation and *Voilà!* — your model is ready to run.

## Getting started


If you are new to AIBECS, head over to the [documentation](https://briochemc.github.io/AIBECS.jl/stable/) and look for the tutorials.
(You can also click on the orange "Documentation" badge above.)

## Concept

This package was developed to exploit linear-algebra tools and algorithms in Julia to efficiently simulate marine tracers.
AIBECS represents global biogeochemical cycles with a discretized system of nonlinear partial differential equations that takes the generic form

(∂/∂𝑡 + 𝓣)*x* = *G*(*x*)

where *x* represents the model state variables, i.e., the marine tracer(s) concentration.
For a single tracer, *x* can be interpreted as the 3D field of its concentration.
In AIBECS, *x* is represented as a column vector.

The operator 𝓣 is a spatial differential operator that represents the transport of tracers.
For example, for a single tracer transported by the ocean circulation,

𝓣 = ∇ ⋅ (***u*** + **K**∇)

represents the effects of advection and eddy-diffusion.
(***u*** is the 3D vector of the marine currents and **K** is a 3×3 diffusivity matrix.)
Thus, 𝓣 *acts* on *x* such that 𝓣*x* is the flux divergence of that tracer.
In AIBECS, 𝓣 is represented by matrices.

Lastly, *G*(*x*) represents the local sources minus sinks of each tracer.
In AIBECS, *G*(*x*) is represented by functions of the tracer(s).

To simulate tracers using the AIBECS, you just need to define the transport operators 𝓣 and the net sources and sinks *G*.
That's pretty much the whole concept!

## References

If you use this package, please cite it.
And if you use data with these package (like the ocean circulation from the OCIM) please also cite them.
The references under bibtex format are available in the [CITATION.bib](./CITATION.bib) file.

Also, if you want to do research using the AIBECS, and you think I could help, do not hesitate to contact me directly (contacts on my [website](www.bpasquier.com)), I would be happy to contribute and collaborate!

<img src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png" alt="NSF" title="NSF_logo" align="right" height="50"/>

The authors acknowledge funding from the Department of Energy grant DE-SC0016539 and from the National Science Foundation grant 1658380.
