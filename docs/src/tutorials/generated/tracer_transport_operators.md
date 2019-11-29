# How the ocean circulation is represented

!!! note
    This tutorial is also available as a Jupyter notebook:
    [`tracer_transport_operators.ipynb`](@__NBVIEWER_ROOT_URL__/tutorials/generated/tracer_transport_operators.ipynb)

## Transport operator

To model marine biogeochemical tracers on a global scale we need to be able to account for their movement within the 3D ocean.
We do this with a tracer **transport operator**, generically denoted by $\mathcal{T}$.
When the operator acts on a tracer field, it produces the divergence of the tracer.
In other words, for a tracer with concentration $x(\boldsymbol{r})$ at location $\boldsymbol{r}$, its divergence at $\boldsymbol{r}$ is $(\mathcal{T} x)(\boldsymbol{r})$.

In the case of the ocean circulation, i.e., the currents and eddies that transport marine tracers floating around in sea water, the transport operator can be represented by

$$\mathcal{T} = \nabla \cdot \left[ \boldsymbol{u}(\boldsymbol{r}) - \mathbf{K}(\boldsymbol{r}) \cdot \nabla \right]$$

where the vector $\boldsymbol{u}$ represents the mean marine current velocity at location $\boldsymbol{r}$.
Thus, $\boldsymbol{u}(\boldsymbol{r})$ is a 3D vector aligned with the ocean currents and whose amplitude, in m/s, gives the speed of the moving sear water at location $\boldsymbol{r}$.
The 3×3 matrix $\mathbf{K}(\boldsymbol{r})$ is the eddy diffusivity matrix.
It represents the effective mixing effect of unresolved eddies, i.e., these turbulent vortices that are too small relative to the model grid to be captured.
Most of the mixing in the ocean occurs along isopycnals, an effect that can be represented by an adequate choice for $\mathbf{K}$.

## Discretization

In order to represent the ocean circulation on a computer one needs to **discretize** the 3D ocean into a discrete grid, i.e., a grid with a finite number of boxes.
One can go a long way towards understanding what a tracer transport operator is by playing with a simple model with only a few boxes, which is the goal of this tutorial.

The simple box model we consider is embeded in a 2×2×2 "shoebox".
It has 5 *wet* boxes and 3 *dry* boxes, as illustrated below:

```@raw html
<img src="https://user-images.githubusercontent.com/4486578/58314610-3b130b80-7e53-11e9-9fe8-9527cdcca2d0.png" width =800>
```

The circulation consists of
- a meridional overturning circulation flowing in a cycle through boxes 1 → 2 → 6 → 5 → 1 (shown in the "meridional section 1" panel above)
- a zonal current in a reentrant cycling through boxes 1 → 3 → 1 (shown in the "layer 1" panel above)
- vertical mixing representing deep convection between boxes 2 ↔ 6 (not shown)

## Vectorization

The 3D tracer field, $x(\boldsymbol{r})$, is **vectorized**, in the sense that the concentrations in each box are rearranged into a column vector.

```@raw html
<img src="https://user-images.githubusercontent.com/4486578/61757212-fe3ba480-ae02-11e9-8d17-d72866eaafb5.gif" width =800>
```

The continuous transport operator $\mathcal{T}$ can then be represented by a matrix, denoted by $\mathbf{T}$, and sometimes called the *transport matrix*.
It turns out that in most cases, this matrix is sparse.

!!! note
    A sparse matrix behaves the same way as a regular matrix.
    The only difference is that in a sparse matrix the majority of the entries are zeros.
    These zeros are not stored explicitly to save computer memory making it possible to deal with fairly high resolution ocean models.

Mathematically, the discretization and vetorization convert an expression with partial derivatives into a matrix vector product.
In summary, for the ocean circulation, we do the following conversion

$$\mathcal{T} x = \nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right] x \longrightarrow \mathbf{T} \, \boldsymbol{x}$$

(We ommitted the $\boldsymbol{r}$ dependency in equations for brevity.)

## Using an ocean circulation in AIBECS

### Loading the grid and the circulation

We will now simulate the transport of a marine tracer controlled by the toy circulation we just saw.
Like for any Julia package, we start by telling Julia that we are going to *use* the package.

```julia
using AIBECS
```

The circulation we just saw is available in AIBECS.
We load it via

```julia
grd, T = Primeau_2x2x2.load() ;
```

which loads 2 objects, `grd`, the model grid, and `T`, the transport matrix.
(Note that in the line above we suppressed the output with `;`, like in MATLAB.)

`grd` is an `OceanGrid` object, of size 2×2×2.

```julia
grd
```

It is a special type for ocean grids, defined in the [OceanGrids](https://github.com/briochemc/OceanGrids.jl) package, on which AIBECS depends.

There are many ways to look inside the grid, one of which is to iterate over it.
For example

```julia
[println(box) for box in grd] ;
```

shows some details about all the boxes of the model, one at a time.
The `grd` object also contains other information about the grid, like the 3D depths of the boxes:

```julia
grd.depth_3D
a
```

or the 3D latitudes:

```julia
grd.lat_3D
```

where you can check that northwards = downwards in the array.
You may notice these come with units!
This helps ensure that degrees of latitude are not confused with meters for example.

!!! note
    Julia depends on and comes with [Unitful](https://github.com/PainterQubits/Unitful.jl), a package for using units, which AIBECS uses.

Inside of `grd`, there is information of which boxes of the grid are "wet" (or "dry").
`wet3D` is a 3D array representing the 3D ocean with `true` for "wet" boxes and `false` for "dry" boxes.
Let's have a look at its contents:

```julia
grd.wet3D
```

It's a 2×2×2 `BitArray`, i.e., an array of bit elements (the `true` and `false` entries).
You can check that it matches our "shoebox" model.
(Well, except for the orientation, for which northwards in the box model is downwards in the array.)

We can find all the wet boxes simply via

```julia
findall(grd.wet3D)
```

These are the 3D indices of the wet boxes, `(i,j,k)`, called the "cartesian" indices.
If you want the "linear" indices, i.e., the numbers as shown in the image of the shoebox model, you can transform `wet3D` into a vector first, via

```julia
iwet = findall(vec(grd.wet3D))
```

We can also check that we indeed have the expected number of wet boxes via

```julia
nb = length(iwet)
```

Now let's have a look at the transport matrix `T`.

```julia
T
```

It's a sparse matrix (of the `SparseMatrixCSC` type in Julia), which only stores non-zero entries.
We can display the full matrix via

```julia
Matrix(T)
```

to check its structure out.

## Radiocarbon

Radiocarbon, ¹⁴C, is produced by cosmic rays in the lower stratosphere and upper troposphere.
It quickly reacts with oxygen to produce ¹⁴CO₂, which is then mixed throughout the troposphere and enters the ocean through air–sea gas exchange.
Because the halflife of radiocarbon is only 5730 years a significant amount of decay can occur before the dissolved inorganic radiocarbon (DI¹⁴C) can mix uniformally throughout the ocean.
As such the ¹⁴C serves as a tracer label for water that was recently in contact with the atmosphere.

### Tracer Equation

Mathematically, the ¹⁴C tracer concentration, denoted $R$ (for Radiocarbon), satisfies the following tracer equation:

$$\frac{\partial R}{\partial t} + \nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right] R = \Lambda(R_\mathsf{atm} - R) - R / \tau,$$

where $\Lambda(R_\mathsf{atm} - R)$ represents the air–sea exchanges and $R / \tau$ the radioactive decay rate.
($\tau$ is the radioactive decay timescale.)

The discretized tracer is thus given by

$$\frac{\partial \boldsymbol{R}}{\partial t} + \mathbf{T} \, \boldsymbol{R} = \mathbf{\Lambda}(R_\mathsf{atm} - \boldsymbol{R}) - \boldsymbol{R} / \tau.$$

### Translation to AIBECS Code

We will perform an idealized radiocarbon simulation in our model and use the ocean circulation defined earlier using AIBECS.
In this model we prescribe the atmospheric concentration, $R_\mathsf{atm}$, to be simply equal to 1.
(We do not specify its unit or its specific value because it is not important for determining the age of a water parcel — only the decay rate does.)

To use AIBECS, one must put the equations into the generic form of

$$\frac{\partial \boldsymbol{x}}{\partial t} + \mathbf{T}(\boldsymbol{p}) \, \boldsymbol{x} = \boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p}),$$

where $\boldsymbol{x}$ is the state vector, $\boldsymbol{p}$ is the vector of model parameters, $\mathbf{T}(\boldsymbol{p})$ is the transport operator, and $\boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p})$ is the local sources minus sinks.

In our radiocarbon-model context, with $\boldsymbol{x} = \boldsymbol{R}$, we have that

$$\boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p}) = \mathbf{\Lambda}(R_\mathsf{atm} - \boldsymbol{x}) - \boldsymbol{x} / \tau.$$

Hence, we must create `T(p)` and `G(x,p)` to give AIBECS the means to simulate the tracer distribution and/or its evolution in time.

#### Sources and Sinks

The local sources and sinks are thus simply given by

```julia
function G(x, p)
    τ, Ratm = p.τ, p.Ratm
    return Λ(Ratm .- x, p) - x / τ
end
```

where `τ` is the decay rate timescale and `Ratm` is the atmospheric concentration of radiocarbon.

We must define the air–sea exchange rate, `Λ(x,p)`, which requires us to define which boxes are located at the surface first.
This is done, e.g., via

```julia
surface_boxes = grd.depth_3D[iwet] .== grd.depth[1]
```

The air–sea exchange rate is then given by

```julia
function Λ(x, p)
    λ, h = p.λ, p.h
    return λ / h * surface_boxes .* x
end
```

where `λ` is the piston velocity and `h` is the height of the top layer of the model grid.

#### Parameters

For the air–sea gas exchange, we use a constant piston velocity $\lambda$ of 50m / 10years, which will happen in the top layer, of height given by, well, the height of the top layer.
And for the radioactive decay we use a timescale $\tau$ of 5730/log(2) years.
We define these as parameters using the dedicated API from the AIBECS:

```julia
t = empty_parameter_table()                   # initialize an empty table of parameters
add_parameter!(t, :τ, 5730u"yr"/log(2)) # radioactive decay e-folding timescale
add_parameter!(t, :λ, 50u"m" / 10u"yr") # piston velocity
add_parameter!(t, :h, grd.δdepth[1])    # height of top layer
add_parameter!(t, :Ratm, 1.0u"mol/m^3") # atmospheric concentration
t
```

shows the parameters that you just created.

We now generate a new object to contain all these parameters via

```julia
initialize_Parameters_type(t, "C14_shoebox_parameters") # creates the type for parameters
p = C14_shoebox_parameters()                            # creates the parameters object
```

#### Generate the state function and its Jacobian

The last step for the setup is for AIBECS to create $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = \boldsymbol{G}(\boldsymbol{x}, \boldsymbol{p}) - \mathbf{T}(\boldsymbol{p}) \, \boldsymbol{x}$, which defines the rate of change of the state, $\boldsymbol{x}$.
This is done via

```julia
F, ∇ₓF = state_function_and_Jacobian(p -> T, G) # generates the state function (and its Jacobian!)
x = zeros(nb)
F(x,p)
```

Note here that AIBECS has automatically created `∇ₓF`, i.e., $\nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p})$, which is the Jacobian of the system.
This Jacobian will be useful in the simulations below.

That's it!
Your model is entirely setup and ready to be used for simulations!

## Time stepping

One way to see how the tracer evolves with time is to time step it.
Here we will use the [Crank-Nicolson](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method) scheme embedded in AIBECS.

A single step is performed via, e.g.,

```julia
δt = ustrip(1.0u"yr" |> u"s")
AIBECS.crank_nicolson_step(x, p, δt, F, ∇ₓF)
```

We can write a function to run all the time steps and save them into a `x_hist` object, via

```julia
function time_steps(x₀, Δt, n, F, ∇ₓF)
    x_hist = [x₀]
    δt = Δt / n
    for i in 1:n
        push!(x_hist, AIBECS.crank_nicolson_step(last(x_hist), p, δt, F, ∇ₓF))
    end
    return reduce(hcat, x_hist), 0:δt:Δt
end
```

!!! note
    We store all the history of `x` in `x_hist`.
    At each step, the `push!` function adds a new `x` to `x_hist`.
    Technically, `x_hist` is a vector of vectors, so at the end, we horizontally concatenate it via `reduce(hcat, x_hist)` to rearrange it into a 2D array.

Now let's simulate the evolution of radiocarbon for 7500 years, starting from a concentration of 1 everywhere, via

```julia
Δt = ustrip(7500u"yr" |> u"s") # 7500 years in seconds
x₀ = ones(5)             #
x_hist, t_hist = time_steps(x₀, Δt, 1000, F, ∇ₓF) # runs the simulation
```

This should take a few seconds to run.

!!! note
    AIBECS also has other time-stepping methods available, including Euler-froward (fast but unstable), Euler-backwards (like Crank-Nicolson, slow but stable), and Crank-Nicolson-Leapfrog (Fast-ish and more stable than Euler forward).

Once it's done, we can plot the evolution of radiocarbon through time via

```julia
using Plots
C14age_hist = -log.(x_hist) * ustrip(p.τ * u"s" |> u"yr")
plt = plot(t_hist * ustrip(1u"s" |> u"yr"), C14age_hist'; label="box " .* string.(iwet))
xlabel!(plt, "simulation time (years)")
ylabel!(plt, "Radiocarbon age (years)")
title!(plt, "Radiocarbon age vs simulation time (Crank-Nicolson)")
```

The box model took more than 4000 years to spin up to equilibrium.
For a box model that's no big deal because it is not computationally expensive to run the model, but for a big circulation model waiting for the model to spinup is painfully long.
We therefore want a better way to find the equilibrium solution.

## Solving Directly for the Steady State

With the AIBECS, you can create the steady state problem and solve it in just 2 commands:

```julia
prob = SteadyStateProblem(F, ∇ₓF, x, p)
x_final = solve(prob, CTKAlg()).u
```

Converting radiocarbon into years gives the following values

```julia
C14age_final = -log.(x_final) * p.τ * u"s" .|> u"yr"
println.("box ", iwet, ": ", C14age_final);
```

These are exactly the limit that the Radiocarbon age reaches after about 4000 years of simulation!
So what happened there?

Well, AIBECS used a version of [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method) to solve for $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0$.
Simply put, Newton's method iterates on the state via the recursion relation

$$\boldsymbol{x}_{n+1} = \boldsymbol{x}_n - \nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p})^{-1} \boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}).$$

Here, our radiocarbon model only requires a single iteration.
This is because the state function, $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p})$, is linear.

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

