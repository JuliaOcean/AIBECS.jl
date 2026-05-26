# AGENTS.md

AIBECS lets you build and solve ocean biogeochemistry models as steady states of

```
∂x/∂t + T(p) * x = G(x, p)
```

where `T(p)` is a sparse transport operator (typically advection + diffusion but can be gravitational settling) and `G(x, p)` is the net local sources and sinks. See [ReadMe.md](ReadMe.md) for the full pitch.

This file has two parts: one for *using* AIBECS, one for *contributing to* it.

---

## 1. Using AIBECS

### The recipe

1. **Load a circulation.** `grd, T = OCIM2.load()` — also available: `OCIM0`, `OCIM1`, `OCIM2_48L`, `OCCA`, `TwoBoxModel`, `Primeau_2x2x2`, `Archer_etal_2000`, `Haine_and_Hall_2025`.
2. **Define `G(x, p)`** — net local sources/sinks as a function of the tracer vector and parameters.
3. **Define a parameter type** `<: AbstractParameters{U}` (a `struct` listing your parameter fields).
4. **Assemble** `F = AIBECSFunction(T, G)` and `prob = SteadyStateProblem(F, x_init, p)`.
5. **Solve** `sol = solve(prob, CTKAlg())`.

### Minimal example (ideal age)

```julia
using AIBECS, JLD2

grd, T = OCIM2.load()
z      = depthvec(grd)

struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end

G(x, p) = @. 1 - x / p.τ * (z ≤ p.z₀)

p      = IdealAgeParameters(1.0, z[1])
F      = AIBECSFunction(T, G)
x_init = zeros(sum(iswet(grd)))
prob   = SteadyStateProblem(F, x_init, p)
age    = solve(prob, CTKAlg())
```

First run will prompt to download the OCIM2 dataset; set `ENV["DATADEPS_ALWAYS_ACCEPT"] = true` to skip the prompt.

### Where to look next

- **Tutorials (start here):** [docs/lit/tutorials/](docs/lit/tutorials/) — `1_ideal_age.jl` is simplest; `3_Pmodel.jl` is a coupled multi-tracer example.
- **How-tos (task-oriented):** [docs/lit/howtos/](docs/lit/howtos/) — parameters, plotting, optimisation, sinking particles, nonlinear solvers.
- **Explanation (concepts):** [docs/src/explanation/](docs/src/explanation/) — transport operators, solvers, datasets.
- **API reference:** [docs/src/reference/functions.md](docs/src/reference/functions.md).
- **Rendered docs:** <https://briochemc.github.io/AIBECS.jl>.

---

## 2. Contributing to AIBECS

### Layout

- [src/AIBECS.jl](src/AIBECS.jl) — module entry point and exports.
- [src/multiTracer.jl](src/multiTracer.jl) — `AIBECSFunction` assembly for single- and multi-tracer models.
- [src/overload_solve.jl](src/overload_solve.jl) — `solve` dispatch for `SteadyStateProblem` + `CTKAlg`.
- [src/CTKsolvers.jl](src/CTKsolvers.jl) — Newton–Chord–Shamanskii solver with line search.
- [src/Parameters.jl](src/Parameters.jl) — `AbstractParameters` machinery (units, bounds, priors).
- [src/CirculationGeneration.jl](src/CirculationGeneration.jl) — `T_advection`, `T_diffusion` helpers.
- Circulation modules in [src/](src/): `OCIM0`, `OCIM1`, `OCIM2`, `OCIM2_48L`, `OCCA`, `TwoBoxModel`, `Primeau_2x2x2`, `Archer_etal_2000`, `Haine_and_Hall_2025`.
- Optional features live in [ext/](ext/) as package extensions (weak deps: `NonlinearSolve`, `NCDatasets`, `JLD2`, `MAT`, `Distributions`, `Bijectors`, `RecipesBase`, `Shapefile`, `Interpolations`).

### Tests

- Entry point: [test/runtests.jl](test/runtests.jl).
- Run locally: `julia --project -e 'using Pkg; Pkg.test()'`.
- Per-feature files in [test/](test/): `setup.jl`, `particles.jl`, `bgc_functions.jl`, `solvers.jl`, `derivatives.jl`, `cost_functions.jl`, …
- CI skips `OCIM2_48L` (memory/time); first-time runs need `DATADEPS_ALWAYS_ACCEPT=true` to auto-fetch circulation data.
- Workflow: [.github/workflows/Test.yml](.github/workflows/Test.yml) — Ubuntu/macOS/Windows × Julia LTS & current.

### Docs build

- Pipeline: Literate.jl → Documenter → DocumenterVitepress, orchestrated by [docs/make.jl](docs/make.jl).
- Sources in [docs/lit/](docs/lit/) (executable `.jl`) compile to Markdown under [docs/src/](docs/src/).
- Fast iteration: `DRAFT=true julia --project=docs docs/make.jl` skips evaluating code blocks.

### Conventions

- Julia ≥ 1.10 (see [Project.toml](Project.toml)).
- Formatting: [Runic.jl](https://github.com/fredrikekre/Runic.jl) — format with `runic -i .` (or check only with `runic --check --diff .`) before committing. Editor integrations are listed in the Runic README.
- No `CONTRIBUTING.md` — for anything not covered above, match surrounding code style.
- Key ecosystem dependencies: `OceanGrids`, `SciMLBase`, `LinearSolve`, `NonlinearSolveBase`, `ForwardDiff`, `Unitful`, `SparseArrays`.
