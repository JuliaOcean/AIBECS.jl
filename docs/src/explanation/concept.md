# [Concept](@id concept)

The AIBECS (pronounced like the cool [ibex](https://en.wikipedia.org/wiki/Ibex) if you have a french accent) is a new software written in [Julia](https://julialang.org) to easily create some marine biogeochemistry models in just a few commands.


AIBECS is not just a single model.
It's a **system** that allows you to create a global steady-state biogeochemistry model with just a few simple commands.
Basically, you just need to tell AIBECS which ocean circulation to use (from simple toy models of just a few boxes to more complicated global models of the circulation), what elements you want to model/track, and how they are converted into other (the net local sources and sinks of your model).
Once the model is set up, chose some parameter values and you can run simulations.

AIBECS relies on many tools from linear algebra to run simulations and perform optimizations really fast.
AIBECS-generated models are described by a state function, denoted $\boldsymbol{F}$, which defines how the concentrations of elements in the ocean evolve with time.
In mathematical terms, this translates to a system of nonlinear differential equations with the generic form 

$$\frac{\partial \boldsymbol{x}}{\partial t} = \boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}),$$
where $\boldsymbol{x}$ is the state of the model (i.e., the concentrations of the tracers), and $\boldsymbol{p}$ are model parameters.
With AIBECS, you can efficiently find the equilibrium of the system (AKA the steady-state).
That is when the time-derivative is $0$, so that

$$\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0,$$

and $\boldsymbol{x}$ does not change with time.
Instead of simulating the evolution of $\boldsymbol{x}$ with time and waiting for the system to reach equilibrium — like most biogeochemistry models do — AIBECS uses linear algebra techniques, like Newton's method in multiple dimensions, or Krylov spaces, to implicitly solve for the steady-state solution, hence the "algebraic" and "implicit" names.
This makes AIBECS much faster than the competition!


