```@raw html
<img src="https://user-images.githubusercontent.com/4486578/60554111-8fc27400-9d79-11e9-9ca7-6d78ee89ea70.png" alt="logo" title="AIBECS_logo" align="right" width="200"/>
```

# AIBECS.jl

**Algebraic Implicit Biogeochemistry Element-Cycling System**

(Work in Progress)

## Introduction

AIBECS may be pronounced like the cool [ibex](https://en.wikipedia.org/wiki/Ibex) (if you have a french accent).
AIBECS is a new software written in [Julia](https://julialang.org) to easily create some marine biogeochmistry models in just a few commands.

!!! note
    AIBECS is developed primarily by Benoît Pasquier with the help of François Primeau and J. Keith Moore from the Department of Earth System Science at the University of California, Irvine, USA.
    This software is in active development, so if you have any suggestions or feature requests, do not hesitate to start an issue directly on the [AIBECS GitHub repository](https://github.com/briochemc/AIBECS.jl)!

AIBECS is not just a single model.
It's a **system** that allows you to create a global steady-state biogeochmistry model with just a few simple commands.
Basically, you just need to tell AIBECS
- (i) which ocean circulation to use (from simple toy models of just a few boxes to more complicated global models of the circulation),
- (ii) what elements you want to model/track and 
- (iii) how each tracer gets converted into other tracers.
Once the model is set up, you can run simulations.

AIBECS relies on many tools from linear algebra to run simulations and perform optimizations really fast.
AIBECS-generated models are described by a state function, denoted $\boldsymbol{F}$, which defines how the concnetrations of elements in the ocean evolve with time.
In mathematical terms, this translates to a system of nonlinear differential equations with the generic form 

$$\frac{\partial \boldsymbol{x}}{\partial t} = \boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}),$$

where $\boldsymbol{x}$ is the state of the model (i.e., the concentrations of the tracers), and $\boldsymbol{p}$ are model parameters.
Here, we are interested in the equilibrium of the system (AKA the steady-state).
That is when the time-derivative is $0$, so that

$$\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0,$$

and $\boldsymbol{x}$ does not change with time.
Instead of simulating the evolution of $\boldsymbol{x}$ with time and waiting for the system to reach equilibrium — like most biogeochemistry models do — AIBECS uses linear algebra techniques, like Newton's method in multiple dimensions, or Krylov spaces, to implicitly solve for the steady-state solution, hence the "algebraic" and "implicit" names.
This makes AIBECS much faster than the competition!

If you want to try AIBECS, head over to the prerequisites page to install the packages (this should take you a few minutes), and then open up one of the notebooks!
