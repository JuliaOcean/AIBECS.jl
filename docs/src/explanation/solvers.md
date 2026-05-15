# [Steady-state solvers](@id solvers)

Finding a steady state of $\boldsymbol{F}(\boldsymbol{x}, \boldsymbol{p}) = 0$
runs an outer **nonlinear solver** that calls an inner **linear solver** at
every Jacobian-update / Newton-step iteration. AIBECS lets you pick each
component independently:

- The outer nonlinear solver can be any algorithm from
  [NonlinearSolve.jl](https://docs.sciml.ai/NonlinearSolve/stable/)
  (`NewtonRaphson`, `TrustRegion`, `LimitedMemoryBroyden`, …), *or*
  AIBECS's own [`CTKAlg`](@ref AIBECS.CTKAlg) — a bespoke
  Newton–Chord–Shamanskii method adapted from C. T. Kelley's MATLAB
  implementations ([Kelley, 2003](https://doi.org/10.1137/1.9780898718898))
  and specialised to the AIBECS state Jacobian.
- The inner linear solver can be any factorisation from
  [LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/) — usually a
  sparse direct factorisation (`UMFPACKFactorization`,
  `KLUFactorization`, `MKLPardisoFactorize`).

The [NonlinearSolve how-to](@ref nonlinearsolve) walks through the call
pattern for the SciML side. For `CTKAlg` the call is just
`solve(prob, CTKAlg())`.

## Recommended pairing: `CTKAlg() + UMFPACKFactorization()`

Calling `CTKAlg()` with no arguments yields
`CTKAlg(linsolve = UMFPACKFactorization())` — the AIBECS-recommended
default. Two reasons it's the recommendation:

1. **CTKAlg uses a Shamanskii update.** The Newton iteration only
   refreshes the Jacobian when the residual stops shrinking fast enough
   (`r > rSham₀ = 0.5`). On AIBECS systems where the steady-state
   Jacobian changes slowly between iterations, one factorisation
   typically covers several Newton steps — generally faster in our
   experience than the "refresh every iteration" baseline that
   `NewtonRaphson` uses. The tradeoff is robustness: when the Jacobian
   does drift, Shamanskii's stale factors can take an extra Newton step
   to recover. Armijo line search is the safety net.
2. **`UMFPACKFactorization` reuses its factors.** `CTKAlg` runs the
   linear solve through a `LinearSolve.LinearCache`, so when Shamanskii
   recycles the Jacobian, UMFPACK's LU factors are reused directly —
   no refactorisation. Pairing CTKAlg with a factorisation that
   supports cache reuse is what makes the Shamanskii win compound.

If `CTKAlg` fails to converge on a pathological problem, the standard
fallback is `NewtonRaphson` (still with UMFPACK) via the NonlinearSolve
route below — fresh Jacobian every iteration, slower but more
forgiving.

## When to use the NonlinearSolve route instead

- You want a different outer algorithm (`TrustRegion`,
  `LimitedMemoryBroyden`, …).
- You want to benchmark several algorithms against `CTKAlg` on your
  own circulation.
- You need NonlinearSolve features that `CTKAlg` does not surface
  (custom termination modes, callbacks, sensitivity hooks).

The cost is the load-time of NonlinearSolve.jl and LinearSolve.jl plus
their transitive deps — neither is a hard AIBECS dependency. The
underlying `SteadyStateProblem` and sparse Jacobian are reused across
both routes; switching is a one-line change at the call site.

## Verified algorithm/factorisation pairs

The combinations below have been checked on the toy `Primeau_2x2x2`
circulation and converge to the same root as `CTKAlg()`.

| Solver call | Match vs `CTKAlg()` | Notes |
|---|---|---|
| `solve(prob, CTKAlg())` | (baseline) | Native AIBECS Newton–Chord–Shamanskii; CTKAlg + UMFPACK default. |
| `solve(nlprob, NewtonRaphson(linsolve = UMFPACKFactorization()))` | exact | Recommended NonlinearSolve fallback. |
| `solve(nlprob, NewtonRaphson(linsolve = KLUFactorization()))` | 1e-10 | KLU may win on very-sparse Jacobians. |
| `solve(nlprob, AIBECS.recommended_nlalg())` | exact | Convenience for `NewtonRaphson + UMFPACK`. |

For the NonlinearSolve rows, `nlprob = AIBECS.nonlinearproblem(prob)` —
see the how-to for the full pattern. The helper attaches the sparse
Jacobian buffer type so LinearSolve's UMFPACK / KLU factorisations
allocate sparse buffers instead of erroring on a dense default.

## Latest benchmark results

Auto-generated from `benchmark/solvers.jl`; the timings below are wall
clock for one steady-state solve on each circulation/tracer pair, on
the maintainer's laptop. The small-tier suite covers `Primeau_2x2x2`,
`OCCA`, and `OCIM0` — enough for relative comparisons across the
available algorithm/factorisation pairs without needing a workstation.

```@eval
import Markdown
raw = read(joinpath(@__DIR__, "latest_timings_small.md"), String)
# Strip the autogen HTML comment, drop the "# Tier: small" metadata
# header (the surrounding section already covers it), then shift all
# remaining headings down one level so they nest under the
# "Latest benchmark results" section.
raw = replace(raw, r"<!--[^\n]*-->\n?" => "")
raw = replace(raw, r"^# Tier:[^\n]*\n+"m => "")
raw = replace(raw, r"^(#+) "m => s"#\1 ")
Markdown.parse(raw)
```

## Krylov solvers

Iterative Krylov backends from LinearSolve (e.g. `KrylovJL_GMRES`,
`KrylovJL_BICGSTAB`) are intentionally not in the supported matrix
above: AIBECS does not ship a preconditioner suited to the
steady-state advection–diffusion–reaction system, and unpreconditioned
Krylov on these operators is uncompetitive with the sparse direct
factorisations. Adding a preconditioner is a future direction.
