# Custom benchmark dashboard — follow-up plan

## Goal

A web dashboard with one chart per `(circulation, tracers)` cell, with one
line per algorithm. e.g. an `OCCA / idealage` chart shows CTKAlg,
NewtonRaphson + UMFPACK, NewtonRaphson + KLU, … as separate lines over
commit history. The use case is spotting when an external solver
suddenly improves or regresses relative to the others — something that
isn't visible across many one-line charts.

## Why we don't already have it

`benchmark-action/github-action-benchmark` renders one chart per unique
bench `name` and has no native multi-series mode. With names like
`"OCCA / idealage / CTKAlg"`, each `(circ, tracer, alg)` triple lands on
its own chart on the stock dashboard.

## Architecture (Oceananigans-style)

Keep `benchmark-action/github-action-benchmark` as the *persistence
pipeline*: it appends a record per commit to `data.js`, handles the
`gh-pages` push, runs the alert threshold, and posts PR comparison
comments. Replace only the *rendering layer* with a custom `index.html`
on the `gh-pages` branch that parses the slash-separated bench names and
draws Chart.js canvases the way we want.

References:

- CliMA/Oceananigans.jl `publish-benchmarks.yml` — workflow that
  transforms benchmark JSON into slash-separated names and hands them
  to github-action-benchmark.
- NumericalEarth/OceananigansBenchmarks `index.html` (on `gh-pages`) —
  the custom Chart.js renderer that reads `data.js` and groups multiple
  bench names onto one chart by parsing the name hierarchy.

Their setup is ~1,700 lines because they have parameter sweeps with
chips, NSYS kernel overlays, animated styling, and a 4-level hierarchy.
Ours needs a 3-level hierarchy (`circ/tracer/alg`) and ~150–250 lines.

## Implementation steps

1. **Bench naming**. In `bench_name(r)` (`benchmark/solvers.jl`),
   switch to slash-only with no surrounding spaces:
   `"OCCA/idealage/CTKAlg"`. The slash becomes the structural separator
   the dashboard JS splits on; the spaces in the current names are only
   there for the MD table (which doesn't need them either now that
   circulation/tracer have moved into section headings).

2. **Workflow**. Revert the per-tier `benchmark-data-dir-path` and
   `name` suffixes in `.github/workflows/Benchmarks.yml`. Tiers control
   what gets *run*, not how the dashboard is *organised* — they
   shouldn't appear in chart paths. Keep a single
   `benchmark-data-dir-path: .` so all tiers' data lands in one
   `data.js`.

3. **Bootstrap `gh-pages`**. Commit a custom `index.html` that:
   - Loads `data.js`.
   - For each `entry.benches[]` splits `name` on `/` into
     `(circ, tracer, alg)`.
   - Builds a hierarchy `(circ, tracer) → alg → [{commit, date, value}]`.
   - Renders one `<canvas>` per `(circ, tracer)` with one Chart.js
     dataset per `alg`, X = commit/date, Y = wall time.
   - Chart.js via CDN
     (`<script src="https://cdn.jsdelivr.net/npm/chart.js@4">`); no
     build step, no npm.

4. **Confirm action doesn't clobber `index.html`**. The action only
   auto-creates `index.html` when absent — so once the custom one is
   committed, subsequent runs should only update `data.js`. Worth
   verifying empirically on the first push; if it does overwrite, add
   a five-line workflow step that restores the custom `index.html`
   after `benchmark-action`.

## When

Defer until both:

- The algorithm matrix, tracer setups, and circulation tiers have
  stabilised. Iterating on those right now would force re-shooting the
  dashboard schema.
- There are ≥10 `main`-branch commits of history. The time-series view
  is the only thing the dashboard adds over the per-tier
  `benchmark/results/latest_timings_${tier}.md` tables, and on a
  near-empty timeline that view is uninformative anyway.

Until then the MD tables cover the side-by-side view and the stock
github-action-benchmark dashboard covers regression alerts.
