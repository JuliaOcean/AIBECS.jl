window.BENCHMARK_DATA = {
  "lastUpdate": 1778340660513,
  "repoUrl": "https://github.com/JuliaOcean/AIBECS.jl",
  "entries": {
    "AIBECS Solver Benchmarks (small)": [
      {
        "commit": {
          "author": {
            "email": "4486578+briochemc@users.noreply.github.com",
            "name": "Benoît Pasquier",
            "username": "briochemc"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b8bc1cf785978170590b97694d37cb3d272afedd",
          "message": "feat: NonlinearSolve.jl + LinearSolve.jl as weak deps (#122)\n\n* feat: NonlinearSolve.jl + LinearSolve.jl as weak deps\n\nAdd an opt-in path to the SciML nonlinear-solver ecosystem alongside\nthe existing CTKAlg native solver. NonlinearSolve, LinearSolve,\nADTypes, DifferentiationInterface, SparseConnectivityTracer, and\nSparseMatrixColorings are wired as weak deps via a new\nAIBECSNonlinearSolveExt extension that exposes:\n\n- AIBECS.nonlinearproblem(prob) — converts a SteadyStateProblem to a\n  NonlinearProblem with sparse jac_prototype attached, sidestepping\n  NonlinearSolve's AutoSpecialize iip-promotion which mishandles\n  AIBECS's 3-arg ODE-style f(u, p, t=0).\n- AIBECS.recommended_nlalg() — NewtonRaphson + UMFPACKFactorization\n  with a sparse-AD fallback (AutoSparse + ForwardDiff via DI +\n  TracerSparsityDetector + GreedyColoringAlgorithm).\n- AIBECS.default_sparse_ad() — exposed for users to override.\n\nCTKAlg, the existing tutorials, and prior SteadyStateProblem usage\nremain unchanged. Tests now run six solver combinations on every\ntoy circulation (CTKAlg + five NonlinearSolve combinations).\n\nAdds a tiered benchmark harness under benchmark/ (toy/small/medium/\nlarge) covering ideal age, radiocarbon, and PO4-POP coupled tracer\nsetups, plus a Benchmarks.yml workflow that publishes to gh-pages\nvia benchmark-action/github-action-benchmark and posts regression\ncomments on PRs. New docs under docs/lit/howtos/nonlinearsolve.jl\nand docs/src/explanation/solvers.md document the integration.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* refactor(NonlinearSolve): lean on documented SciML API, drop AD scaffolding\n\nThe original PR shipped a manual SteadyStateProblem to NonlinearProblem shim\nthat pre-dated the current SciMLBase.NonlinearProblem(::SteadyStateProblem)\nauto-conversion, plus a sparse-AD fallback (ADTypes, DifferentiationInterface,\nSparseConnectivityTracer, SparseMatrixColorings) that is redundant given\nAIBECS always supplies an analytical Jacobian.\n\nThis commit:\n- Shrinks ext/AIBECSNonlinearSolveExt to a thin wrapper around the SciMLBase\n  conversion plus a jac_prototype attachment so LinearSolve allocates a\n  sparse buffer (UMFPACK/KLU otherwise error on the dense default).\n- Drops ADTypes / DifferentiationInterface / SparseConnectivityTracer /\n  SparseMatrixColorings from weakdeps, compat, extras, [extensions], targets.\n- Drops TrustRegion / LevenbergMarquardt rows from the test and benchmark\n  matrices: on AIBECS-style near-linear systems they stall via\n  ReturnCode.ShrinkThresholdExceeded (catastrophic cancellation in the\n  trust-region acceptance ratio near the root) at residual ~1e-11 while\n  NewtonRaphson collapses to ~1e-19. Not a different root, just lower\n  precision; not worth supporting.\n- Adds NonlinearSolve and LinearSolve to docs/Project.toml so the new\n  how-to @example blocks evaluate.\n- Switches benchmark/tracer_setups.jl to `import AIBECS: @units, units`\n  matching the convention used in src/Parameters.jl, test/parameters.jl,\n  and the docs/lit tutorials (`using` cannot extend `units`).\n- Drops the stale [sources] block from benchmark/Project.toml that pinned\n  AIBECS to a different machine dev path. The CI workflow already calls\n  Pkg.develop(path=pwd()), and locally the same pattern works.\n- Adds MUMPS_jll, Pardiso, ParU_jll, Sparspak to benchmark deps in\n  preparation for a wider linear-solver matrix (probe work in progress;\n  ParU_jll currently ships without OpenMP so it runs serial).\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* benchmark: align tracer setups with tutorials, drop redundant ustrip wrapper\n\nThe bundled circulation loaders in ext/AIBECSJLD2Ext.jl already strip\nunits at load time, so the `_strip` wrapper in benchmark/circulations.jl\nwas a no-op and is removed.\n\nThe Pmodel benchmark setup had drifted from docs/lit/tutorials/3_Pmodel.jl\n(hand-rolled S₀/S′ via PFDO instead of `transportoperator`, ztop==0 mask\ninstead of (z ≤ z₀), missing -POP/τ_geo term, hard-coded literal in\nsms_POP, divergent parameter struct field names). Rewritten to mirror\nthe tutorial verbatim: same equations, same parameter values, same\nPmodelParameters struct (with @initial_value defaults), same x_init.\n\nIdeal age and radiocarbon setups already matched their tutorials and\nare unchanged.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* style: split multi-package using lines + apply Runic formatting\n\n- Each `using A, B, C` is split into one `using` per line so that adding\n  or removing a single dependency is a single-line diff. Selective\n  imports (`using A: x, y`) are left grouped since they pull from one\n  package.\n- `julia -m Runic -i` applied across src/, ext/, test/, benchmark/,\n  docs/make.jl, and docs/lit/.\n- Move each `export ModuleName` from the bottom of its module file into\n  src/AIBECS.jl alongside the corresponding `include`. Without this\n  Runic indents the module body (see fredrikekre/Runic.jl#5) because\n  there is content after the closing `end`.\n\nNo behavioural changes; full test suite still passes.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* benchmark(fix): load JLD2 so OCIM*/OCCA loaders work in CI\n\nOCIM0.load and friends are stubs in src/ that error until\nAIBECSJLD2Ext is activated. Adding `using JLD2` to benchmark/solvers.jl\nfixes the `using JLD2` requirement check that the Benchmarks.yml\nworkflow was hitting.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* benchmark: fix NLsolveJL stall, restructure outputs, add dashboard plan\n\nThe shared `abstol = eps(Float64)` was an unreachable convergence target\nfor the wrapper algorithms (NLsolveJL, SIAMFANLEquationsJL). Native\nNewtonRaphson overrides convergence via the RelNormSafeBest\ntermination_condition; the wrappers reject that kwarg and use abstol as\nthe literal |F| tolerance. At n≈10⁵ matvec roundoff alone exceeds eps,\nso NLsolveJL silently spun to maxiters=1000 — the \"hang\" the prior\nanalysis pinned on the linsolve. Split a wrappers-only kwargs bundle\n(abstol=1e-8, maxiters=50, show_trace).\n\nDrop NLSolversJL: its NonlinearSolve extension allocates a dense N×N\nJacobian (~57 GB at OCCA scale). Keep the LinearSolve.jl/UMFPACK\nadapter on NLsolveJL so the wrapper exercises nonlinear logic against\nthe same factorization the native row uses.\n\nTier-suffixed output files (latest_timings_\\${tier}.md /\nbench_results_\\${tier}.json) so successive runs of different tiers\naccumulate locally instead of clobbering. MD groups by (circulation,\ntracers), integer seconds. CI workflow gains per-tier artifact names,\ndashboard data dirs, and docs sync paths; BENCH_TIER promoted to\nworkflow-level env so all steps see it.\n\nAlign tutorial/howto z₀ with the harness (z[1] not 30.0).\n\nbenchmark/DASHBOARD_PLAN.md sketches the Oceananigans-style multi-line\nper-(circ, tracer) chart we deferred — github-action-benchmark for data\npersistence, custom Chart.js index.html for rendering.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* benchmark(plan): note unconditional PR-comment posting of latest timings\n\nThe current PR-comparison step only fires on regression past\n`alert-threshold`, so on a fresh repo with no gh-pages baseline the\nbenchmark feels invisible during review. Capture this as follow-up\nin DASHBOARD_PLAN.md so we remember to wire up an always-post step\n(e.g. via a sticky comment) independently of the Chart.js work.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* benchmark: cap iterations at 20 + always post timings on PRs\n\nCap `maxiters = 20` on every solver row (CTKAlg's `maxItNewton = 20`,\nnative NewtonRaphson, and the SIAMFANLEquationsJL/NLsolveJL wrappers).\nAt AIBECS's tracer setups anything that hasn't converged in 20 Newton\nsteps is not going to, so this trades a non-convergence flag for\nbounded wall time and bounded log size — last CI run hit the wrappers'\ndefault 1000-iter cap and the runner timed out after 2h47m before\nuploading any logs.\n\nAlso add a `marocchino/sticky-pull-request-comment` step that posts\n`latest_timings_${tier}.md` on every PR push, keyed by tier, so the\ntimings table is visible during review even when the comparison step\nfinds nothing past `alert-threshold`. Trim the now-implemented Aside\nout of `DASHBOARD_PLAN.md` and replace it with a one-paragraph\npointer to the workflow step.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* benchmark: scale-aware non-dim wrapper + per-component convergence test\n\nAdd `benchmark/dimensional_scaling.jl` implementing the Dennis-Schnabel\n`(D_x, D_f)` diagonal scaling (DOI 10.1137/1.9781611971200, Ch. 7) as a\nlocal helper: `nondimensionalize(prob; scale_u, scale_F)` returns a\n`ScaledProblem` whose `.scaled` field is the rescaled\n`SteadyStateProblem` (state and residual O(1) per component). The\nJacobian transform uses `lmul!`/`rmul!` against `Diagonal` for in-place\nsparse row+column scaling, preserving `SparseMatrixCSC{Float64,Int64}`\nso UMFPACK/KLU/ParU keep their fast paths. `default_scales` provides a\nKINSOL-style fallback; `pmodel_scales` gives the physics-based recipe\nfor the coupled DIP/POP system from steady-state mass balance.\n\nRewire `solvers.jl` to non-dimensionalise once per `(circulation,\ntracer)` pair and feed every algorithm the scaled problem. Replace the\n`reltol = 1/τstop ≈ 3e-14` bopping pathology with a scalar\n`NL_ABSTOL = 1e-6` against the Inf-norm of the scaled residual,\naligning every solver's internal norm to Inf where configurable\n(CTKAlg via `nrm`, NewtonRaphson via `internalnorm`); SIAMFANL is\nalready Inf, NLsolveJL is stuck at L2 because the SciML extension\nhardcodes it (caveat surfaced in the per-solve diagnostic). Each row\nprints both the universal `max-rel-drift = max_i |F_i/scale_F_i|` and\nthat solver's own internal stopping quantity, evaluated on the final\niterate.\n\nReplace the residual MD column with `Max rel drift` (dimensionless,\ncomparable across `(circ, tracer)` setups); JSON `extra` field follows.\n\nAlso start `po4pop` from `vcat(DIP_geo·ones(nb), zeros(nb))` instead\nof `DIP_geo·ones(2nb)` — the previous initial guess over-seeded POP\nby ~46x its steady-state scale (≈ DIP_geo · τ_POP/τ_DIP), driving a\nhuge first Newton step. Zero POP is closer to the true steady state.\n\nNote `DASHBOARD_PLAN.md` aside on the column meaning change.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* Fix OCCA mass conservation, add AD + conservation tests, refactor\n\nOCCA non-convergence (visible in last benchmark CI: all six solvers stalled\non `OCCA / po4pop`) traced to a non-mass-conserving raw transport operator\nthat wasn't being corrected on load. Fixed alongside several supporting\npieces from earlier sessions; bundling them here.\n\nOCCA grid:\n- Refresh OCCA dataset to a new Zenodo release (URL + MD5 bumped).\n- Apply default `conservemass = true` correction in `OCCA.load`: subtract\n  `Diagonal((Tᵀ v) ./ v)` so `vᵀT = 0` by construction. Pass\n  `conservemass = false` to recover the raw matrix; in that case tracer\n  budgets drift unless the geological restoring is modeled explicitly.\n- Doc/warning block describing the choice.\n\nDataDeps staleness:\n- New `src/datadeps_cache.jl` with `_invalidate_stale_cache(name, filename,\n  expected_md5)` helper that walks the DataDeps load path and removes a\n  cached directory whose MD5 no longer matches what's registered. Closes\n  the gap where `@datadep_str` short-circuits without rechecking, so\n  bumping a registered URL/MD5 in source has no effect on stale clones.\n- Hooked into every JLD2-backed `*.load()` (`OCIM0`, `OCIM1`, `OCIM2`,\n  `OCCA`) at the top of the loader.\n\nAD tests:\n- New `test/sparse_jacobian.jl` validates AIBECS's bespoke\n  `local_jacobian`-assembled `prob.f.jac` against a reference Jacobian\n  computed via DifferentiationInterface + SparseConnectivityTracer +\n  ForwardDiff. Run for all three benchmark problems (idealage,\n  radiocarbon, po4pop) on every circulation in `test_setup_only` and on\n  every toy circulation in `test_everything`.\n- `runtests.jl` now `include`s `benchmark/problems.jl` so the test suite\n  shares the same problem builders as the benchmark harness.\n\nConservation tests:\n- `test/setup.jl` adds two checks per circulation: `T · 1 = 0`\n  (divergence-free) and `vᵀ T = 0` (mass-conserving). OCCA's\n  `conservemass=true` enforces the latter at the cost of perturbing\n  column sums, so `T · 1 = 0` is `@test_broken` for OCCA only.\n\nBenchmark refactor:\n- Extract `IdealAgeParameters`, `RadiocarbonParameters`, `PmodelParameters`,\n  the matching `build_*_problem` builders, and `pmodel_scales` from\n  `benchmark/solvers.jl` into `benchmark/problems.jl`, included near the\n  top alongside `dimensional_scaling.jl`. Same file is reused by\n  `test/sparse_jacobian.jl` so the model setup is single-sourced.\n- `solvers.jl` shrinks by ~110 lines.\n\nMisc:\n- `Project.toml` cleanup.\n- Formatting and typo fixes across `src/{ETOPO, GroundWaters, OCCA, OCIM*,\n  Primeau_2x2x2, aeolian_sources}.jl`, `ext/AIBECSOCIM2_48LExt.jl`.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* benchmark: warmup + single timed solve, drop SIAMFANL/ParU, prune deps\n\nRestructure the per-row measurement protocol and prune the algorithm\nmatrix to keep CI wall time honest.\n\nMeasurement protocol:\n- Each `(circ, tracer, alg)` row now does an untimed warmup solve on a\n  fresh problem (pays first-call JIT compile), rebuilds, then runs a\n  single `@timed solve` on the new fresh problem. The previous \"1st +\n  2nd timed solve\" pattern showed essentially no difference between\n  passes in practice (caches the algorithms keep on `prob` between\n  calls don't reduce wall time meaningfully here), so the second\n  timing is dropped along with its column in the MD table and the\n  `1st_solve` field in the JSON `extra`.\n- Output collapses to a single `Time (s)` column.\n\nAlgorithm pruning:\n- Drop `NewtonRaphson + ParU`: hangs intermittently inside METIS\n  NodeND on OCIM0-class sparsity (`ParU_C_Analyze →\n  umf_l_cholmod → cholmod_l_metis → recursive\n  MlevelNestedDissection`). Reproduced again on the latest small-tier\n  CI run.\n- Drop `SIAMFANLEquationsJL`: ~10× slower first-call compile than\n  NewtonRaphson, and no per-iter trace exposed by the SciML wrapper\n  (only an end-of-solve `printerr` for failures). Removed for clarity\n  rather than a correctness reason.\n- `ALG_MATRIX` is now CTKAlg, NewtonRaphson + UMFPACK,\n  NewtonRaphson + KLU, NLsolveJL + Newton/UMFPACK. Top-of-file comment\n  block documents why each removed solver is excluded.\n- Drop the corresponding `import SIAMFANLEquations` and `import\n  ParU_jll` from `solvers.jl` and `format_diagnostic`'s SIAMFANL\n  branch.\n\nPer-tracer timescale:\n- Print `τ = ‖x‖_v / ‖F‖_v` per tracer block on every converged row,\n  where `‖y‖_v = √(yᵀ Diagonal(v) y / Σv)` is the volume-weighted RMS\n  (volumes from `OceanGrids.volumevec(grd)`). Computed in *original*\n  units so the result is in seconds; printed via `t * u\"s\" |> u\"yr\"`\n  pipe so the `yr` label is genuinely meaningful. Used only for\n  display; the convergence test still uses the unweighted Inf-norm of\n  the scaled residual (`max-rel-drift ≤ NL_ABSTOL = 1e-6`).\n- New `tracer_names` field on each `TRACER_SETUPS` entry labels the\n  per-tracer τ in the output (`age` for idealage, `R` for radiocarbon,\n  `DIP, POP` for po4pop). The vector's length determines how the\n  state vector is sliced into per-tracer blocks.\n- New `τ (yr)` column in the MD table.\n\nProject.toml:\n- Drop unused deps: `BenchmarkTools`, `ForwardDiff`, `MUMPS_jll`,\n  `NCDatasets`, `Pardiso`, `SIAMFANLEquations`, `ParU_jll`, `Sparspak`,\n  `Statistics`. None are imported by `benchmark/*.jl`.\n- Add `Dates` (was used via stdlib transient resolution but not\n  declared).\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* test: add @info banners, isolate sparse Jacobian to OCCA\n\nPer-file @info lines make CI logs interpretable when a runner is killed\nmid-test (the symptom on Linux Julia 1/lts in the prior run).\n\nRestrict the SCT+ForwardDiff sparse-Jacobian test to OCCA — OCIM1's\nreference Jacobian build was too heavy for GitHub-hosted Linux runners.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* test: fix Circulation binding in OCCA sparse-Jacobian testset\n\nWas assigning the Symbol `:OCCA` instead of the `OCCA` module, so\n`Circulation.load()` failed with \"type Symbol has no field `load`\".\nThe for-loop sites use `$C` interpolation to splice the bare name in;\nhere we just write `OCCA` directly.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\n\n* CTKAlg: pluggable LinearSolve, NonlinearSolveBase termination, Pardiso bench rows\n\n- CTKAlg{L} carries a linsolve field forwarded to LinearSolve. The\n  Newton-Chord-Shamanskii inner loop now mutates a LinearSolve.LinearCache\n  (linsol.A= / linsol.b=) instead of calling factorize and backslash, so the\n  Shamanskii lazy refresh maps directly onto LinearSolve's isfresh toggle.\n  AgedJacobianFactors and its backslash overload are gone - age is a plain Int.\n- Convergence delegates to NonlinearSolveBase termination modes (default\n  AbsNormSafeBestTerminationMode(maximum-of-abs), matching NonlinearSolve's\n  :regular default. abstol/reltol flow through; pass any\n  AbstractNonlinearTerminationMode to swap criteria. Replaces the old\n  nrm/tau-stop kwargs.\n- Promotes LinearSolve from weakdep to direct dep; adds NonlinearSolveBase.\n- Benchmark gains CTKAlg + UMFPACK, CTKAlg + KLU, and 4 conditional MKL\n  Pardiso rows (NewtonRaphson x 2, CTKAlg x 2) gated on\n  Pardiso.mkl_is_available() with nprocs = Sys.CPU_THREADS (override via\n  BENCH_NPROCS). CTK now consumes the same AbsNormSafeBest termination\n  as NewtonRaphson, so format_diagnostic collapses the CTK branch.\n- Verified locally on macOS: Pkg.test 95/95 pass; small-tier benchmark shows\n  CTKAlg and NewtonRaphson produce identical residuals at every linsolve at\n  OCIM0 scale (n = 382338). Pardiso rows skip silently on macOS as designed\n  and will activate on the ubuntu-latest runner.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>\nEOF\n)\n\n---------\n\nCo-authored-by: Claude Opus 4.7 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-05-09T13:48:08+10:00",
          "tree_id": "b1b83312ddf66b7634859d3c6328b0fe217b8023",
          "url": "https://github.com/JuliaOcean/AIBECS.jl/commit/b8bc1cf785978170590b97694d37cb3d272afedd"
        },
        "date": 1778299961351,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "OCCA / idealage / CTKAlg + UMFPACK",
            "value": 0.737538789,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.067e-11; τ_yr=age=8.29e+14"
          },
          {
            "name": "OCCA / idealage / CTKAlg + KLU",
            "value": 4.211708776,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=6.821e-11; τ_yr=age=5.28e+14"
          },
          {
            "name": "OCCA / idealage / NewtonRaphson + UMFPACK",
            "value": 0.993414255,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.067e-11; τ_yr=age=8.29e+14"
          },
          {
            "name": "OCCA / idealage / NewtonRaphson + KLU",
            "value": 3.915365976,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=6.821e-11; τ_yr=age=5.28e+14"
          },
          {
            "name": "OCCA / idealage / NLsolveJL + Newton/UMFPACK",
            "value": 0.993608789,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.067e-11; τ_yr=age=8.29e+14"
          },
          {
            "name": "OCCA / idealage / NewtonRaphson + MKLPardisoFactorize",
            "value": 0.680734872,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.276e-12; τ_yr=age=2.59e+15"
          },
          {
            "name": "OCCA / idealage / NewtonRaphson + MKLPardisoIterate",
            "value": 0.70957973,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.276e-12; τ_yr=age=2.59e+15"
          },
          {
            "name": "OCCA / idealage / CTKAlg + MKLPardisoFactorize",
            "value": 0.671211678,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.276e-12; τ_yr=age=2.59e+15"
          },
          {
            "name": "OCCA / idealage / CTKAlg + MKLPardisoIterate",
            "value": 0.672794709,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.276e-12; τ_yr=age=2.59e+15"
          },
          {
            "name": "OCCA / radiocarbon / CTKAlg + UMFPACK",
            "value": 0.693977799,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.098e-08; τ_yr=R=2.24e+13"
          },
          {
            "name": "OCCA / radiocarbon / CTKAlg + KLU",
            "value": 3.913366106,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.357e-08; τ_yr=R=1.71e+13"
          },
          {
            "name": "OCCA / radiocarbon / NewtonRaphson + UMFPACK",
            "value": 0.978343212,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.098e-08; τ_yr=R=2.24e+13"
          },
          {
            "name": "OCCA / radiocarbon / NewtonRaphson + KLU",
            "value": 4.224109239,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.357e-08; τ_yr=R=1.71e+13"
          },
          {
            "name": "OCCA / radiocarbon / NLsolveJL + Newton/UMFPACK",
            "value": 1.050750622,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.098e-08; τ_yr=R=2.24e+13"
          },
          {
            "name": "OCCA / radiocarbon / NewtonRaphson + MKLPardisoFactorize",
            "value": 0.744321437,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=2.235e-08; τ_yr=R=3.89e+13"
          },
          {
            "name": "OCCA / radiocarbon / NewtonRaphson + MKLPardisoIterate",
            "value": 0.748787637,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=2.235e-08; τ_yr=R=3.89e+13"
          },
          {
            "name": "OCCA / radiocarbon / CTKAlg + MKLPardisoFactorize",
            "value": 1.01519105,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=2.235e-08; τ_yr=R=3.89e+13"
          },
          {
            "name": "OCCA / radiocarbon / CTKAlg + MKLPardisoIterate",
            "value": 0.718467774,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=2.235e-08; τ_yr=R=3.89e+13"
          },
          {
            "name": "OCCA / po4pop / CTKAlg + UMFPACK",
            "value": 1.808610298,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.498e-07; τ_yr=DIP=5e+09,POP=2.5e+06"
          },
          {
            "name": "OCCA / po4pop / CTKAlg + KLU",
            "value": 4.661535467,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.498e-07; τ_yr=DIP=5e+09,POP=2.5e+06"
          },
          {
            "name": "OCCA / po4pop / NewtonRaphson + UMFPACK",
            "value": 3.063613672,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.811e-10; τ_yr=DIP=3.84e+12,POP=1.93e+09"
          },
          {
            "name": "OCCA / po4pop / NewtonRaphson + KLU",
            "value": 11.235197845,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.811e-10; τ_yr=DIP=3.84e+12,POP=1.93e+09"
          },
          {
            "name": "OCCA / po4pop / NLsolveJL + Newton/UMFPACK",
            "value": 3.425793271,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.811e-10; τ_yr=DIP=3.84e+12,POP=1.93e+09"
          },
          {
            "name": "OCCA / po4pop / NewtonRaphson + MKLPardisoFactorize",
            "value": 4.093394629,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.811e-10; τ_yr=DIP=3.84e+12,POP=1.93e+09"
          },
          {
            "name": "OCCA / po4pop / NewtonRaphson + MKLPardisoIterate",
            "value": 4.19069364,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.811e-10; τ_yr=DIP=3.84e+12,POP=1.93e+09"
          },
          {
            "name": "OCCA / po4pop / CTKAlg + MKLPardisoFactorize",
            "value": 1.980583108,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.498e-07; τ_yr=DIP=5e+09,POP=2.5e+06"
          },
          {
            "name": "OCCA / po4pop / CTKAlg + MKLPardisoIterate",
            "value": 1.993367816,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.498e-07; τ_yr=DIP=5e+09,POP=2.5e+06"
          },
          {
            "name": "OCIM0 / idealage / CTKAlg + UMFPACK",
            "value": 4.205739047,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=5.675e-10; τ_yr=age=1.29e+14"
          },
          {
            "name": "OCIM0 / idealage / CTKAlg + KLU",
            "value": 47.704789581,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.060e-09; τ_yr=age=8.37e+13"
          },
          {
            "name": "OCIM0 / idealage / NewtonRaphson + UMFPACK",
            "value": 4.144900465,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=5.675e-10; τ_yr=age=1.29e+14"
          },
          {
            "name": "OCIM0 / idealage / NewtonRaphson + KLU",
            "value": 47.006195112,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.060e-09; τ_yr=age=8.37e+13"
          },
          {
            "name": "OCIM0 / idealage / NLsolveJL + Newton/UMFPACK",
            "value": 4.594531746,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=5.675e-10; τ_yr=age=1.29e+14"
          },
          {
            "name": "OCIM0 / idealage / NewtonRaphson + MKLPardisoFactorize",
            "value": 8.255122032,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.507e-11; τ_yr=age=1.23e+15"
          },
          {
            "name": "OCIM0 / idealage / NewtonRaphson + MKLPardisoIterate",
            "value": 8.482962153,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.078e-11; τ_yr=age=1.26e+15"
          },
          {
            "name": "OCIM0 / idealage / CTKAlg + MKLPardisoFactorize",
            "value": 4.194725329,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.507e-11; τ_yr=age=1.23e+15"
          },
          {
            "name": "OCIM0 / idealage / CTKAlg + MKLPardisoIterate",
            "value": 4.269606293,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.078e-11; τ_yr=age=1.26e+15"
          },
          {
            "name": "OCIM0 / radiocarbon / CTKAlg + UMFPACK",
            "value": 4.019984523,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.134e-09; τ_yr=R=8.65e+13"
          },
          {
            "name": "OCIM0 / radiocarbon / CTKAlg + KLU",
            "value": 47.214995273,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.136e-08; τ_yr=R=5.49e+13"
          },
          {
            "name": "OCIM0 / radiocarbon / NewtonRaphson + UMFPACK",
            "value": 4.485674342,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.134e-09; τ_yr=R=8.65e+13"
          },
          {
            "name": "OCIM0 / radiocarbon / NewtonRaphson + KLU",
            "value": 46.395819738,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.136e-08; τ_yr=R=5.49e+13"
          },
          {
            "name": "OCIM0 / radiocarbon / NLsolveJL + Newton/UMFPACK",
            "value": 4.296897548,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.134e-09; τ_yr=R=8.65e+13"
          },
          {
            "name": "OCIM0 / radiocarbon / NewtonRaphson + MKLPardisoFactorize",
            "value": 8.573281448,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.062e-09; τ_yr=R=4.28e+14"
          },
          {
            "name": "OCIM0 / radiocarbon / NewtonRaphson + MKLPardisoIterate",
            "value": 8.917176148,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.022e-09; τ_yr=R=4.3e+14"
          },
          {
            "name": "OCIM0 / radiocarbon / CTKAlg + MKLPardisoFactorize",
            "value": 4.582990789,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.062e-09; τ_yr=R=4.28e+14"
          },
          {
            "name": "OCIM0 / radiocarbon / CTKAlg + MKLPardisoIterate",
            "value": 4.65662048,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.022e-09; τ_yr=R=4.29e+14"
          },
          {
            "name": "OCIM0 / po4pop / CTKAlg + UMFPACK",
            "value": 6.343493977,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=5.733e-09; τ_yr=DIP=3.68e+10,POP=2.9e+07"
          },
          {
            "name": "OCIM0 / po4pop / CTKAlg + KLU",
            "value": 47.063950545,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=5.733e-09; τ_yr=DIP=3.68e+10,POP=2.9e+07"
          },
          {
            "name": "OCIM0 / po4pop / NewtonRaphson + UMFPACK",
            "value": 9.163952993,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=1.017e-08; τ_yr=DIP=2.48e+10,POP=1.96e+07"
          },
          {
            "name": "OCIM0 / po4pop / NewtonRaphson + KLU",
            "value": 89.531531442,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=1.017e-08; τ_yr=DIP=2.48e+10,POP=1.96e+07"
          },
          {
            "name": "OCIM0 / po4pop / NLsolveJL + Newton/UMFPACK",
            "value": 11.561464741,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=1.017e-08; τ_yr=DIP=2.48e+10,POP=1.96e+07"
          },
          {
            "name": "OCIM0 / po4pop / NewtonRaphson + MKLPardisoFactorize",
            "value": 23.804133253,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=1.010e-07; τ_yr=DIP=3.9e+09,POP=2.71e+10"
          },
          {
            "name": "OCIM0 / po4pop / NewtonRaphson + MKLPardisoIterate",
            "value": 17.068750448,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=2.419e-07; τ_yr=DIP=1.85e+09,POP=1.96e+07"
          },
          {
            "name": "OCIM0 / po4pop / CTKAlg + MKLPardisoFactorize",
            "value": 9.199858826,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=7.898e-07; τ_yr=DIP=5.25e+08,POP=9.43e+06"
          },
          {
            "name": "OCIM0 / po4pop / CTKAlg + MKLPardisoIterate",
            "value": 10.107859535,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=5.733e-09; τ_yr=DIP=3.67e+10,POP=2.9e+07"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "briochemc@gmail.com",
            "name": "Benoit Pasquier",
            "username": "briochemc"
          },
          "committer": {
            "email": "briochemc@gmail.com",
            "name": "Benoit Pasquier",
            "username": "briochemc"
          },
          "distinct": true,
          "id": "c0899e091f630d10e1dfcf1879fa2aef6561c305",
          "message": "docs: switch NonlinearSolve howto to Primeau_2x2x2 to fix CI runtime\n\nThe 3D OCIM0 grid combined with NonlinearSolve's NewtonRaphson + UMFPACK\ntook ~89 minutes per @example block on macOS docs CI (vs ~3 min for\nCTKAlg on the same problem), pushing total docs build time from ~25 min\nto 3-8 hours. Primeau_2x2x2 matches the \"verified on\" scope already\nstated in the solvers explanation page and runs both paths in seconds.\n\nBumps patch to 0.17.1.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-05-10T01:06:37+10:00",
          "tree_id": "ee88afaab06639c627f2b1dbc3e49957e1f308b2",
          "url": "https://github.com/JuliaOcean/AIBECS.jl/commit/c0899e091f630d10e1dfcf1879fa2aef6561c305"
        },
        "date": 1778340653793,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "OCCA / idealage / CTKAlg + UMFPACK",
            "value": 1.14458778,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.067e-11; τ_yr=age=8.29e+14"
          },
          {
            "name": "OCCA / idealage / CTKAlg + KLU",
            "value": 4.703810509,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=6.821e-11; τ_yr=age=5.28e+14"
          },
          {
            "name": "OCCA / idealage / NewtonRaphson + UMFPACK",
            "value": 1.117859316,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.067e-11; τ_yr=age=8.29e+14"
          },
          {
            "name": "OCCA / idealage / NewtonRaphson + KLU",
            "value": 4.993035856,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=6.821e-11; τ_yr=age=5.28e+14"
          },
          {
            "name": "OCCA / idealage / NLsolveJL + Newton/UMFPACK",
            "value": 1.112713842,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.067e-11; τ_yr=age=8.29e+14"
          },
          {
            "name": "OCCA / idealage / NewtonRaphson + MKLPardisoFactorize",
            "value": 0.809736396,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.276e-12; τ_yr=age=2.59e+15"
          },
          {
            "name": "OCCA / idealage / NewtonRaphson + MKLPardisoIterate",
            "value": 1.115685367,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.276e-12; τ_yr=age=2.59e+15"
          },
          {
            "name": "OCCA / idealage / CTKAlg + MKLPardisoFactorize",
            "value": 0.765902681,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.276e-12; τ_yr=age=2.59e+15"
          },
          {
            "name": "OCCA / idealage / CTKAlg + MKLPardisoIterate",
            "value": 1.092177546,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.276e-12; τ_yr=age=2.59e+15"
          },
          {
            "name": "OCCA / radiocarbon / CTKAlg + UMFPACK",
            "value": 1.043072488,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.098e-08; τ_yr=R=2.24e+13"
          },
          {
            "name": "OCCA / radiocarbon / CTKAlg + KLU",
            "value": 4.67006156,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.357e-08; τ_yr=R=1.71e+13"
          },
          {
            "name": "OCCA / radiocarbon / NewtonRaphson + UMFPACK",
            "value": 0.709824001,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.098e-08; τ_yr=R=2.24e+13"
          },
          {
            "name": "OCCA / radiocarbon / NewtonRaphson + KLU",
            "value": 5.013485922,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=7.357e-08; τ_yr=R=1.71e+13"
          },
          {
            "name": "OCCA / radiocarbon / NLsolveJL + Newton/UMFPACK",
            "value": 0.712547158,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=4.098e-08; τ_yr=R=2.24e+13"
          },
          {
            "name": "OCCA / radiocarbon / NewtonRaphson + MKLPardisoFactorize",
            "value": 1.115573861,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=2.235e-08; τ_yr=R=3.89e+13"
          },
          {
            "name": "OCCA / radiocarbon / NewtonRaphson + MKLPardisoIterate",
            "value": 0.81495899,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=2.235e-08; τ_yr=R=3.89e+13"
          },
          {
            "name": "OCCA / radiocarbon / CTKAlg + MKLPardisoFactorize",
            "value": 0.812751714,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=2.235e-08; τ_yr=R=3.89e+13"
          },
          {
            "name": "OCCA / radiocarbon / CTKAlg + MKLPardisoIterate",
            "value": 0.771509776,
            "unit": "s/solve",
            "extra": "n=84661; max_rel_drift=2.235e-08; τ_yr=R=3.89e+13"
          },
          {
            "name": "OCCA / po4pop / CTKAlg + UMFPACK",
            "value": 2.229662276,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.498e-07; τ_yr=DIP=5e+09,POP=2.5e+06"
          },
          {
            "name": "OCCA / po4pop / CTKAlg + KLU",
            "value": 4.711067098,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.498e-07; τ_yr=DIP=5e+09,POP=2.5e+06"
          },
          {
            "name": "OCCA / po4pop / NewtonRaphson + UMFPACK",
            "value": 3.352808151,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.811e-10; τ_yr=DIP=3.84e+12,POP=1.93e+09"
          },
          {
            "name": "OCCA / po4pop / NewtonRaphson + KLU",
            "value": 12.834645845,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.811e-10; τ_yr=DIP=3.84e+12,POP=1.93e+09"
          },
          {
            "name": "OCCA / po4pop / NLsolveJL + Newton/UMFPACK",
            "value": 3.888737644,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.811e-10; τ_yr=DIP=3.84e+12,POP=1.93e+09"
          },
          {
            "name": "OCCA / po4pop / NewtonRaphson + MKLPardisoFactorize",
            "value": 4.966446532,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.811e-10; τ_yr=DIP=3.84e+12,POP=1.93e+09"
          },
          {
            "name": "OCCA / po4pop / NewtonRaphson + MKLPardisoIterate",
            "value": 5.049972261,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.811e-10; τ_yr=DIP=3.84e+12,POP=1.93e+09"
          },
          {
            "name": "OCCA / po4pop / CTKAlg + MKLPardisoFactorize",
            "value": 2.32156835,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.498e-07; τ_yr=DIP=5e+09,POP=2.5e+06"
          },
          {
            "name": "OCCA / po4pop / CTKAlg + MKLPardisoIterate",
            "value": 2.232158477,
            "unit": "s/solve",
            "extra": "n=169322; max_rel_drift=1.498e-07; τ_yr=DIP=5e+09,POP=2.5e+06"
          },
          {
            "name": "OCIM0 / idealage / CTKAlg + UMFPACK",
            "value": 4.401463508,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=5.675e-10; τ_yr=age=1.29e+14"
          },
          {
            "name": "OCIM0 / idealage / CTKAlg + KLU",
            "value": 59.475586813,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.060e-09; τ_yr=age=8.37e+13"
          },
          {
            "name": "OCIM0 / idealage / NewtonRaphson + UMFPACK",
            "value": 3.970082646,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=5.675e-10; τ_yr=age=1.29e+14"
          },
          {
            "name": "OCIM0 / idealage / NewtonRaphson + KLU",
            "value": 58.995545745,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.060e-09; τ_yr=age=8.37e+13"
          },
          {
            "name": "OCIM0 / idealage / NLsolveJL + Newton/UMFPACK",
            "value": 4.113943059,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=5.675e-10; τ_yr=age=1.29e+14"
          },
          {
            "name": "OCIM0 / idealage / NewtonRaphson + MKLPardisoFactorize",
            "value": 9.363189391,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.507e-11; τ_yr=age=1.23e+15"
          },
          {
            "name": "OCIM0 / idealage / NewtonRaphson + MKLPardisoIterate",
            "value": 9.720425576,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.078e-11; τ_yr=age=1.26e+15"
          },
          {
            "name": "OCIM0 / idealage / CTKAlg + MKLPardisoFactorize",
            "value": 4.803431701,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.507e-11; τ_yr=age=1.23e+15"
          },
          {
            "name": "OCIM0 / idealage / CTKAlg + MKLPardisoIterate",
            "value": 5.325718354,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.078e-11; τ_yr=age=1.26e+15"
          },
          {
            "name": "OCIM0 / radiocarbon / CTKAlg + UMFPACK",
            "value": 4.299307258,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.134e-09; τ_yr=R=8.65e+13"
          },
          {
            "name": "OCIM0 / radiocarbon / CTKAlg + KLU",
            "value": 58.937064902,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.136e-08; τ_yr=R=5.49e+13"
          },
          {
            "name": "OCIM0 / radiocarbon / NewtonRaphson + UMFPACK",
            "value": 4.737759472,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.134e-09; τ_yr=R=8.65e+13"
          },
          {
            "name": "OCIM0 / radiocarbon / NewtonRaphson + KLU",
            "value": 58.05687255,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.136e-08; τ_yr=R=5.49e+13"
          },
          {
            "name": "OCIM0 / radiocarbon / NLsolveJL + Newton/UMFPACK",
            "value": 4.717411255,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=3.134e-09; τ_yr=R=8.65e+13"
          },
          {
            "name": "OCIM0 / radiocarbon / NewtonRaphson + MKLPardisoFactorize",
            "value": 9.920166757,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.062e-09; τ_yr=R=4.28e+14"
          },
          {
            "name": "OCIM0 / radiocarbon / NewtonRaphson + MKLPardisoIterate",
            "value": 10.26721098,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.022e-09; τ_yr=R=4.3e+14"
          },
          {
            "name": "OCIM0 / radiocarbon / CTKAlg + MKLPardisoFactorize",
            "value": 4.905959784,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.062e-09; τ_yr=R=4.28e+14"
          },
          {
            "name": "OCIM0 / radiocarbon / CTKAlg + MKLPardisoIterate",
            "value": 5.04597602,
            "unit": "s/solve",
            "extra": "n=191169; max_rel_drift=1.022e-09; τ_yr=R=4.29e+14"
          },
          {
            "name": "OCIM0 / po4pop / CTKAlg + UMFPACK",
            "value": 7.220192779,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=5.733e-09; τ_yr=DIP=3.68e+10,POP=2.9e+07"
          },
          {
            "name": "OCIM0 / po4pop / CTKAlg + KLU",
            "value": 59.576350629,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=5.733e-09; τ_yr=DIP=3.68e+10,POP=2.9e+07"
          },
          {
            "name": "OCIM0 / po4pop / NewtonRaphson + UMFPACK",
            "value": 10.819169513,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=1.017e-08; τ_yr=DIP=2.48e+10,POP=1.96e+07"
          },
          {
            "name": "OCIM0 / po4pop / NewtonRaphson + KLU",
            "value": 111.724492161,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=1.017e-08; τ_yr=DIP=2.48e+10,POP=1.96e+07"
          },
          {
            "name": "OCIM0 / po4pop / NLsolveJL + Newton/UMFPACK",
            "value": 13.112701381,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=1.017e-08; τ_yr=DIP=2.48e+10,POP=1.96e+07"
          },
          {
            "name": "OCIM0 / po4pop / NewtonRaphson + MKLPardisoFactorize",
            "value": 27.48236299,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=1.010e-07; τ_yr=DIP=3.9e+09,POP=2.71e+10"
          },
          {
            "name": "OCIM0 / po4pop / NewtonRaphson + MKLPardisoIterate",
            "value": 19.634436192,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=2.419e-07; τ_yr=DIP=1.85e+09,POP=1.96e+07"
          },
          {
            "name": "OCIM0 / po4pop / CTKAlg + MKLPardisoFactorize",
            "value": 10.644073606,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=7.898e-07; τ_yr=DIP=5.25e+08,POP=9.43e+06"
          },
          {
            "name": "OCIM0 / po4pop / CTKAlg + MKLPardisoIterate",
            "value": 11.406543405,
            "unit": "s/solve",
            "extra": "n=382338; max_rel_drift=5.733e-09; τ_yr=DIP=3.67e+10,POP=2.9e+07"
          }
        ]
      }
    ]
  }
}