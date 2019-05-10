
<img src="https://user-images.githubusercontent.com/4486578/57189839-c422db00-6f56-11e9-9e1d-26c8d9208702.png" alt="logo" title="AIBECS_logo" align="right" height="400"/>




# AIBECS

*The ideal tool for exploring global marine biogeochemical cycles.*


<p>
  <img src="https://img.shields.io/badge/stability-experimental-orange.svg">
</p>
<p>
  <a href="https://doi.org/<DOI_NUMBER>">
    <img src="https://zenodo.org/badge/DOI/<DOI_NUMBER>.svg" alt="DOI">
  </a>
  <a href="https://github.com/briochemc/AIBECS.jl/blob/master/LICENSE">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-yellow.svg">
  </a>
</p>
<p>
  <a href="https://briochemc.github.io/AIBECS.jl/stable/">
    <img src=https://img.shields.io/badge/docs-stable-blue.svg>
  </a>
  <a href="https://briochemc.github.io/AIBECS.jl/latest/">
    <img src=https://img.shields.io/badge/docs-dev-blue.svg>
  </a>
</p>
<p>
  <a href="https://travis-ci.com/briochemc/AIBECS.jl">
    <img alt="Build Status" src="https://travis-ci.com/briochemc/AIBECS.jl.svg?branch=master">
  </a>
  <a href='https://coveralls.io/github/briochemc/AIBECS.jl'>
    <img src='https://coveralls.io/repos/github/briochemc/AIBECS.jl/badge.svg' alt='Coverage Status' />
  </a>
</p>
<p>
  <a href="https://ci.appveyor.com/project/briochemc/AIBECS-jl">
    <img alt="Build Status" src="https://ci.appveyor.com/api/projects/status/prm2xfd6q5pba1om?svg=true">
  </a>
  <a href="https://codecov.io/gh/briochemc/AIBECS.jl">
    <img src="https://codecov.io/gh/briochemc/AIBECS.jl/branch/master/graph/badge.svg" />
  </a>
</p>






**AIBECS** (for **A**lgebraic **I**mplicit **B**iogeochemical **E**lemental **C**ycling **S**ystem, pronounced like the cool [ibex](https://en.wikipedia.org/wiki/Ibex)) is a Julia package that provides ocean biogeochmistry modelers with an easy-to-use interface for creating and running models of the ocean system.

AIBECS is a system because it allows you to chose some biogeochemical tracers, define their interactions, select an ocean circulation and *Voilà!* — your model is ready to run.

## Getting started

We are currently working on the documentation along with some simple Jupyter notebooks to get you started with AIBECS:
- modeling the ideal mean age (amlost ready)
- ventilation tracers (coming soon)
- coupled nutrients (coming soon)


## Motivation

The idea for this package came about in part from the [AWESOME OCIM](https://github.com/hengdiliang/AWESOME-OCIM-v1.1) by [Seth John](https://dornsife.usc.edu/cf/earth/faculty_display.cfm?Person_ID=1063621) and others.
The idea behind the AWESOME OCIM is that modeling simple global steady-state marine biogeochemical tracers should be an easy task.
The AWSEOME OCIM provides a MATLAB GUI to model biogeochemical tracers embedded in a sparse "transport" matrix circulation, AKA the **O**cean **C**irculation **I**nverse **M**odel (OCIM) by [Tim DeVries](https://tdevries.eri.ucsb.edu).
(OCIM matrices and references can be found on Tim's website [here](https://tdevries.eri.ucsb.edu/models-and-data-products/).)
However awesome it is, the AWESOME OCIM lacks some important features, which motivated this package.
A non-exhaustive list of some of the features that we were looking for:
- Use something else than MATLAB, which is not available to everyone.
    [Julia](https://julialang.org) provides the perfect free alternative, with better performance and better syntax.
- Allow for nonlinear systems.
    OCIM users regularly model nonlinear mechanisms, and use Newton-type solvers to run simulations.
    The AIBECS makes it easy to model nonlinear mechanisms by providing solvers under the hood, so that you don't have to worry about them.
    Additionally, having these solvers in an open package that is thoroughly tested greatly reduces the chance of bugs!
- Allow for coupling of tracers, which is fundamental to our understanding of global marine biogeochmical cycles.
    The AIBECS aims to provide the easiest possible interface for you to create multi-tracer models.
    In the tests (and soon in the documentation), should be some examples of multiple and nonlinear tracer model implementations.
- Allow optimizations of model parameters.
    Arguably, using the fast simulations that are afforded by steady-state circulations should standardize objective optimization of model parameters constrained by available observational data (see, e.g., [Pasquier and Holzer, 2017](https://www.biogeosciences.net/14/4125/2017/)).
    With AIBECS, we developed a state-of-the-art autodifferentation tool, the [F-1 method](https://github.com/briochemc/F1Method.jl) (see Pasquier et al., in preparation).
    It was developed specifically for this type of optimizations to run as fast as possible, i.e., it allows you to compute gradient and Hessians of an objective function as fast as if you had gone through the trouble of deriving each second-order derivative by hand! 
    (In the future, monthly circulation matrices, see, e.g., the CYCLOCIM project, should be available from AIBECS) 
- Use other circulations than the output of the OCIM.
    AIBECS aims to provide a simple API for you to load any transport matrix for the ocean circulation, as long as the matrix creator made it publicly available.
- Plotting publication-quality figures.
    Every modeler has reinvented the wheel when it comes to plotting.
    But it should not be this way.
    Eventually, AIBECS will provide users with a plotting interface that can be used to directly produce flawless vectorized figures for your publications.

We emphasize that this package is under active development, so that not all the features advertized above are implemented.
(Plotting publication-quality figures will likely be a feature that takes time, considering the current state of plotting in Julia!)


## The Maths

In AIBECS, global biogeochemical cycles are represented by discretized nonlinear partial differential equations that take the generic form

```julia
∂x/∂t = F(x,p)
```

where `x` is a column vector of the model state variables (i.e., the tracers) and `p` is a vector of model parameters.
(For now, AIBECS only handles steady models, for which `F` does not depend on time.)

This package was developed for models to exploit techniques from linear algebra.
A typical example is if the model is linear (or affine to be specific), i.e., if

```julia
F(x,p,t) = A * x + b
```

then the model's steady state can be computed in a single use of "backslash", via `s = A \ -b`.
(That's what the AWESOME OCIM does.)

However, AIBECS also works for nonlinear steady-state problems, i.e., where `F(x,p)` is nonlinear, covering a much larger applications!
In this case, AIBECS can use a state-of-the-art Newton type of solver to find the steady-state solution for you.
(The solver was adapted from the quasi-Newton solver written in MATLAB by C.T. Kelley)

## References

Please cite us of you use this package.
After the first release, their should be a citable DOI produced by Zenodo.
The reference information wil also eventually be available in the CITATION.bib file.
