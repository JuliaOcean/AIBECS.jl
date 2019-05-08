
<img src="https://user-images.githubusercontent.com/4486578/57189839-c422db00-6f56-11e9-9e1d-26c8d9208702.png" alt="logo" title="AIBECS_logo" align="center" height="400"/>


The ideal tool for exploring global marine biogeochemical cycles.


# AIBECS

**AIBECS** (for **A**lgebraic **I**mplicit **B**iogeochemical **E**lemental **C**ycling **S**ystem, pronounced like the cool [ibex](https://en.wikipedia.org/wiki/Ibex)) is a Julia package that provides ocean biogeochmistry modelers with an easy-to-use interface for creating and running models of the ocean system.

AIBECS is a system because it allows you to chose some biogeochemical tracers, define their interactions, select an ocean circulation and *Voilà!* — your model is ready to run.

## Getting started

We have prepared basic Jupyter notebooks to get you started with AIBECS:
- modeling the ideal mean age

## Motivation

Rsearchers spend a significant amount of time developing models or tools to answer scientific questions.
Under the publish-or-perish pressure, they often develop subpar code that either already exists somewhere else, contains mistakes, or simply does not work as they think it does.
A solution is to package reusable code as open-source software to perform those tasks, hopefully allowing researchers to focus on science rather than on technical development.
Said code should be easy to read, fast, open-source, and available to anyone to use for free — for which [Julia](https://julialang.org) is perfectly suited!

The motivation for this particular package came in part from the development of the [AWESOME OCIM](https://github.com/hengdiliang/AWESOME-OCIM-v1.1) by [Seth John](https://dornsife.usc.edu/cf/earth/faculty_display.cfm?Person_ID=1063621) and others.
The idea behind the AWESOME OCIM is that modeling simple global steady-state marine biogeochemical tracers should be an easy task.
The AWSEOME OCIM uses a sparse "transport" matrix to represent the ocean circulation from the output of the **O**cean **C**irculation **I**nverse **M**odel (OCIM) by [Tim DeVries](https://tdevries.eri.ucsb.edu).
(The OCIM matrices and references can be found on Tim's website [here](https://tdevries.eri.ucsb.edu/models-and-data-products/).)
However awesome it is the AWESOME OCIM lacks a few good things, which motivated this package.
A non-exhaustive list of these reasons is listed here:
- it uses MATLAB, which is not available to everyone. 
    Julia provides the perfect free alternative, with better performance and better syntax. 
- it allows only for linear systems. 
    But OCIM users regularly model nonlinear mechanisms, and use Newton-type solvers to run simulations.
    AIBECS uses such solvers under the hood, without you worrying about them.
    Additionally, having these solvers in an open package that is thoroughly tested greatly reduces the chance of bugs.
- it does not allow for coupling of tracers, which is fundamental to our understanding of global marine biogeochmical cycles.
    AIBECS aims to provide the easiest possible interface for you to create multi-tracer models.
- it does not allow you to run optimizations of your model parameters.
    Arguably, using the fast simulations that are afforded by steady-state circulations (or seasonally steady ones like the matrices under development in the CYCLOCIM project) should standardize objective optimization of model parameters constrained by available observational data (see, e.g., [Pasquier and Holzer, 2017](https://www.biogeosciences.net/14/4125/2017/))
    With AIBECS, we developed a state-of-the-art autodifferentation tool, the [F-1 method](https://github.com/briochemc/F1Method.jl) (see Pasquier et al., in preparation).
    It was developed specifically for this type of optimizations to run as fast as possible, i.e., it allows you to compute gradient and Hessians of an objective function as fast as if you had gone through the trouble of deriving each second-order derivative by hand!
- it only uses the OCIM circulation.
    AIBECS aims to provide a single API to load any transport matrix for the ocean circulation, as long as the matrix creator make it available publicly.
- it uses prescribed observational data that may have been preprocessed in ways you would not want it to be.
    Instead, AIBECS aims to allow you to download raw observational data and process it as you wish to constrain your model.
- its plots are not publication-ready.
    Eventually, AIBECS will provide its users with a plotting interface that can be used to directly produce flawless vectorized figures for your publications.

We emphasize that this package is under active development, so that not all the features advertized above are implemented.
(Plotting publication-quality figures will likely be a feature that takes time, considering the current state of plotting in Julia.)


## The Maths

In AIBECS, global biogeochemical cycles are represented by discretized nonlinear partial differential equations that take the generic form

<img src="https://latex.codecogs.com/svg.latex?&space;\frac{\partial&space;\boldsymbol{x}}{\partial&space;t}&space;=&space;\boldsymbol{F}(\boldsymbol{x},&space;\boldsymbol{p},&space;t)" title="Eq1"/>

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\boldsymbol{x}" title="\boldsymbol{x}" /> is a column vector of the model state variables (i.e., the tracers) and <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\boldsymbol{p}" title="p"/> is a vector of model parameters.

This package was developed for models to exploit techniques from linear algebra.
A typical example is if the model is linear, i.e, if

<img src="https://latex.codecogs.com/svg.latex?&space;\boldsymbol{F}(\boldsymbol{x},&space;\boldsymbol{p},&space;t)&space;=&space;\mathbf{A}\,\boldsymbol{x}" title="Eq2"/>

then the model's steady state can be computed in a single matrix inversion.
However, AIBECS also works for nonlinear steady-state problems, i.e., where

<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{F}(\boldsymbol{x},&space;\boldsymbol{p})&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\boldsymbol{F}(\boldsymbol{x},&space;\boldsymbol{p})&space;=&space;0" title="\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0" /></a>

In this case, AIBECS can use a state-of-the-art Newton type of solver to find the steady-state solution for you.

## Features / TODO list

List of features to be implemented.
This is work in progress and serves me as a TODO list as I develop the package.

-  [x] Load OCIM's transport matrix and grid information.
-  [x] Develop a biogeochemical model for tracer(s)—This task inherently must be the user's responsibility. However it may be possible to simplify this part of the workflow.
-  [ ] Steady-state solvers:
    - [x] Newton for cases where `F` is nonlinear and small (solver code adapted from C.T. Kelley (2003))
    - [ ] Newton-Krylov if nonlinear and big.
-  [ ] Time-stepper.
-  [x] Optimize biogeochemical parameters by minimizing some objective function (has been moved out of the package)
-  [ ] Generate linear equivalent model and diagnostics ([Pasquier and Holzer, 2018](https://www.biogeosciences.net/15/7177/2018/)).
-  [ ] Plot stuff.


## References

Please cite us of you use this package.
If available, use the CITATION.bib file.
