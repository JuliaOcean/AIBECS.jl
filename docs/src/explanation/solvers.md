# [Steady-state solvers](@id solvers)

AIBECS exposes two routes to the steady state of a tracer system
$\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0$:

1. **`CTKAlg()`** — AIBECS's native Newton-Krylov solver, ported from
   C. T. Kelley's textbook implementations. It ships inside AIBECS,
   pulls in no extra dependencies, and is the algorithm used by every
   tutorial.
2. **NonlinearSolve.jl + LinearSolve.jl** — an opt-in integration that
   gives access to the broader [SciML](https://sciml.ai) nonlinear
   solver ecosystem. The wiring lives in `ext/AIBECSNonlinearSolveExt.jl`
   and only loads once the user runs `using NonlinearSolve, LinearSolve`.

The [NonlinearSolve how-to](@ref nonlinearsolve) walks through the call
pattern; this page explains *when* you might want one path over the other
and lists the algorithm/factorisation pairs that have been verified.

## When to use which

- **Use `CTKAlg()` if** you want the smallest possible install, you're
  running the tutorials as written, or you need a stable baseline to
  compare against. It is the zero-dependency default.
- **Use NonlinearSolve if** you want to swap in a different LinearSolve
  factorisation (UMFPACK, KLU, …) or benchmark several algorithms
  against `CTKAlg`. Activating the integration costs the load-time of
  NonlinearSolve.jl and LinearSolve.jl plus their transitive deps.

Switching between the two paths is a one-line change at the call site —
the underlying `SteadyStateProblem` and its sparse Jacobian are reused.

## Algorithms

The combinations below have been verified on the toy `Primeau_2x2x2`
circulation and converge to the same root as `CTKAlg()`.

| Solver call | Match vs `CTKAlg()` | Notes |
|---|---|---|
| `solve(prob, CTKAlg())` | (baseline) | Native AIBECS Newton-Krylov; consumes the `SteadyStateProblem`. |
| `solve(nlprob, NewtonRaphson(linsolve = UMFPACKFactorization()))` | exact | Recommended default. |
| `solve(nlprob, NewtonRaphson(linsolve = KLUFactorization()))` | 1e-10 | KLU may be faster on very-sparse Jacobians. |
| `solve(nlprob, AIBECS.recommended_nlalg())` | exact | NewtonRaphson + UMFPACK. |

For the NonlinearSolve rows, `nlprob = AIBECS.nonlinearproblem(prob)` —
see the how-to for the full pattern. The helper attaches the sparse
Jacobian buffer type so LinearSolve's UMFPACK / KLU factorisations
allocate sparse buffers instead of erroring on a dense default.

## Latest benchmark results (auto-updated)

A placeholder until the benchmark workflow lands and starts committing
fresh numbers.

| Circulation | Tracers | Algorithm | n | Time (s) | Allocs (MB) | Residual | Converged |
|---|---|---|---:|---:|---:|---|:---:|
| _Pending first benchmark run_ | | | | | | | |

```@eval
# Future: replace this stub with a Documenter `include` of
# docs/src/explanation/latest_timings.md once the benchmark workflow
# commits it. For now, the table above is a static placeholder.
```

## Recommendations

Recommendations will firm up once benchmark numbers are in. In the
meantime:

- **UMFPACK is the safe sparse default.** Mature, well-tested, and the
  factorisation used by the AIBECS-recommended algorithm.
- **KLU may win on very-sparse Jacobians.** Worth trying if a profile
  shows the linear solve dominating.

## Krylov solvers

Iterative Krylov backends from LinearSolve (e.g. `KrylovJL_GMRES`,
`KrylovJL_BICGSTAB`) are intentionally not in the supported matrix yet:
AIBECS does not currently ship a preconditioner suited to the
steady-state advection–diffusion–reaction system, and unpreconditioned
Krylov on these operators is uncompetitive with the sparse direct
factorisations above. Adding a preconditioner is a future direction.
