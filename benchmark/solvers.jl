# Benchmark harness for AIBECS steady-state solvers.
#
# Run with:
#   julia --project=benchmark benchmark/solvers.jl
#
# Configuration via env vars:
#   BENCH_TIER    ∈ {toy, small, medium, large}            (default: small)
#   BENCH_TRACERS ∈ {idealage, radiocarbon, po4pop, all}   (default: all)
#
# Goal: help users pick a solver + linear-factorisation pair on
# AIBECS-shaped problems. Output is intended for embedding in
# `docs/src/explanation/solvers.md`, not for commit-over-commit perf
# tracking.
#
# Two orthogonal sweeps, both anchored at the same reference algorithm
# `REFAlg = CTKAlg(linsolve = UMFPACKFactorization())` (the AIBECS-recommended
# default) so each sweep's plot reads as "what changes if I vary just this
# axis from the recommended setup?":
#
#   - LINEAR_SWEEP   — vary the inner LinearSolve factorisation (UMFPACK,
#                      KLU, MKL Pardiso) under fixed CTKAlg nonlinear.
#   - NONLINEAR_SWEEP — vary the outer nonlinear algorithm (CTKAlg,
#                      NewtonRaphson, NLsolveJL) under fixed UMFPACK linear.
#
# `REFAlg` is the first row in both sweeps. It is timed once and rendered
# in both output sections so each sweep can be read standalone.
#
# Outputs:
#   - benchmark/results/latest_timings_<tier>.md
#       Two sub-tables per (circ, tracer), one per sweep. Consumed by
#       `docs/make.jl` via the `@eval` block in `solvers.md`.
#   - docs/src/explanation/figures/bench_<sweep>_<circ>_<tracer>_<tier>.png
#       2x1 layout per cell (walltime bar chart on top, iter/Jacobian/
#       linear-solve diagnostics scatter below, x-axis-linked). One PNG
#       per (sweep, circ, tracer) — small/focused figures rather than
#       one mega-plot per sweep.
#
# Scaling:
#   The harness also runs each cell twice — once with `nondimensionalize`
#   applied, once unscaled. The scaled pass is **skipped** when the tracer
#   setup does not supply its own `scales` recipe (`scales = nothing`):
#   the fallback `default_scales` heuristic in
#   `benchmark/dimensional_scaling.jl` is degenerate when `u0 ≈ 0` (it
#   produces `scale_u ≈ eps` everywhere → ũ ~ 1/eps), so scaling those
#   cells says nothing useful about whether scaling helps. Only `po4pop`
#   currently ships physics-based scales and gets both rows.
#
# Solvers intentionally omitted:
# - NLSolversJL: its NonlinearSolve extension doesn't plumb `jac_prototype`
#   through to NLSolvers.jl, which then allocates `J = similar(x, N, N)` dense
#   (~57 GB at n≈84k). Skip until upstream lands.
# - SIAMFANLEquationsJL: very slow first-call compile (~10× slower than
#   NewtonRaphson) and no per-iter trace. Removed for benchmark clarity.
# - NewtonRaphson + ParU: hangs inside METIS NodeND on OCIM0 sparsity.
# - Krylov: AIBECS lacks a preconditioner suitable for the steady-state
#   advection-diffusion-reaction system.

using AIBECS
using JLD2          # activates AIBECSJLD2Ext so OCIM*/OCCA loaders work
using NonlinearSolve
using NonlinearSolveBase
using LinearSolve
using SparseArrays
using Printf
using LinearAlgebra
using SciMLBase
using Dates
using UnPack
using Unitful
using Unitful: m, d, yr, Myr, mmol, μmol
using CairoMakie
import AIBECS: @units, units
import AIBECS: @initial_value, initial_value
import NLsolve               # activates NonlinearSolveNLsolveExt
import Pardiso               # activates LinearSolvePardisoExt; needed for MKLPardiso* constructors

CairoMakie.activate!(type = "png")

include(joinpath(@__DIR__, "dimensional_scaling.jl"))
include(joinpath(@__DIR__, "problems.jl"))

# === Circulation tiers =====================================================
# `toy` runs purely-bundled circulations (no DataDeps, no JLD2 download). Best
# for iterating on the harness itself or sanity-checking solver behaviour.

const TIERS = (
    toy = [
        ("Primeau_2x2x2",     Primeau_2x2x2.load),
        ("TwoBoxModel",       TwoBoxModel.load),
        ("Archer_etal_2000",  Archer_etal_2000.load),
        ("Haine_and_Hall_2025", Haine_and_Hall_2025.load),
    ],
    small = [("OCCA", OCCA.load), ("OCIM0", OCIM0.load)],
    medium = [("OCIM1", OCIM1.load), ("OCIM2", OCIM2.load)],
    large = [("OCIM2_48L", OCIM2_48L.load)],
)

# Tiny helper so progress reaches the terminal / log file as it happens —
# `@info` is line-buffered when the script is piped, which makes the run
# look like it has stalled.
function vprintln(args...)
    println(args...)
    flush(stdout)
    return nothing
end

# === Tracer setups (mirror the corresponding tutorials) ====================
# `IdealAgeParameters`, `RadiocarbonParameters`, `PmodelParameters` and the
# matching `build_*_problem` builders + `pmodel_scales` live in
# `benchmark/problems.jl`. The same file is reused by `test/sparse_jacobian.jl`.
# `scales = nothing` ⇒ the heuristic `default_scales` would be used, but
# because `u0 = zeros` for the single-tracer setups that heuristic is
# degenerate (see top-of-file note); the harness only runs the *scaled*
# pass when `scales !== nothing`.
const TRACER_SETUPS = (
    idealage    = (label = "idealage",    build = build_idealage_problem,    scales = nothing,       tracer_names = ["age"]),
    radiocarbon = (label = "radiocarbon", build = build_radiocarbon_problem, scales = nothing,       tracer_names = ["R"]),
    po4pop      = (label = "po4pop",      build = build_pmodel_problem,      scales = pmodel_scales, tracer_names = ["DIP", "POP"]),
)

# === Algorithm sweeps ======================================================
# `family` classifies the wrapper kind for the stats-extraction routine.
# CTKAlg surfaces stats via `sol.stats::SciMLBase.NLStats` (populated in
# `src/overload_solve.jl`); NewtonRaphson surfaces them via the same field;
# NLsolveJL wraps the original NLsolve result on `sol.original`.

# NLsolveJL's default `linsolve = (x, A, b) -> copyto!(x, A \ b)` falls
# through to Julia's generic sparse solve. Routing the per-step solve through
# LinearSolve.jl / UMFPACK matches NewtonRaphson + UMFPACK's linear-algebra
# path so the wrapper benchmark exercises NLsolve's nonlinear logic against
# the same factorization the native row uses.
nlsolve_linsolve_umfpack(x, A, b) =
    (copyto!(x, solve(LinearProblem(A, b), UMFPACKFactorization()).u); x)

# MKL Pardiso rows are platform-conditional. `Pardiso.mkl_is_available()`
# returns true on Linux (where MKL_jll ships native binaries) and false on
# macOS, so dev runs on a Mac silently skip those rows. `nprocs` is the
# OMP-thread count for Pardiso; on `ubuntu-latest` runners `Sys.CPU_THREADS`
# picks up the runner's full core count. `BENCH_NPROCS` overrides.
const PARDISO_AVAILABLE = Pardiso.mkl_is_available()
const PARDISO_NPROCS = parse(Int, get(ENV, "BENCH_NPROCS", string(Sys.CPU_THREADS)))
if PARDISO_AVAILABLE
    vprintln("MKL Pardiso available — adding Pardiso rows with nprocs = $PARDISO_NPROCS")
else
    vprintln("MKL Pardiso not available on this platform — skipping Pardiso rows")
end

# Reference algorithm — the AIBECS recommendation. Both sweeps anchor here
# so the cross-sweep comparison is fair.
const REFAlg = (
    label = "CTKAlg + UMFPACK",
    alg = CTKAlg(linsolve = UMFPACKFactorization()),
    family = :ctk, needs_nlprob = false,
)

# Linear-solver sweep — CTKAlg fixed, `linsolve` varied.
const LINEAR_SWEEP = [
    REFAlg,
    (label = "CTKAlg + KLU", alg = CTKAlg(linsolve = KLUFactorization()), family = :ctk, needs_nlprob = false),
    (PARDISO_AVAILABLE ? [
        (label = "CTKAlg + MKLPardisoFactorize", alg = CTKAlg(linsolve = MKLPardisoFactorize(; nprocs = PARDISO_NPROCS)), family = :ctk, needs_nlprob = false),
        (label = "CTKAlg + MKLPardisoIterate",   alg = CTKAlg(linsolve = MKLPardisoIterate(; nprocs = PARDISO_NPROCS)),   family = :ctk, needs_nlprob = false),
    ] : [])...,
]

# Nonlinear-solver sweep — UMFPACK fixed, outer nonlinear algorithm varied.
const NONLINEAR_SWEEP = [
    REFAlg,
    (label = "NewtonRaphson + UMFPACK",    alg = NewtonRaphson(linsolve = UMFPACKFactorization()),                       family = :nr,      needs_nlprob = true),
    (label = "NLsolveJL + Newton/UMFPACK", alg = NLsolveJL(method = :newton, linsolve = nlsolve_linsolve_umfpack),        family = :nlsolve, needs_nlprob = true),
]

const SWEEPS = (
    linear = (
        title = "Linear-solver sweep",
        subtitle = "CTKAlg held fixed; `linsolve` varied",
        cases = LINEAR_SWEEP,
    ),
    nonlinear = (
        title = "Nonlinear-solver sweep",
        subtitle = "UMFPACK held fixed; nonlinear algorithm varied",
        cases = NONLINEAR_SWEEP,
    ),
)

# Deduplicated list of cases actually run end-to-end. `REFAlg` is the
# shared baseline (same `label`), so we dedupe by label to time it once.
const ALL_CASES = let
    cases = NamedTuple[]
    seen = Set{String}()
    for case in vcat((s.cases for s in SWEEPS)...)
        case.label in seen && continue
        push!(seen, case.label)
        push!(cases, case)
    end
    cases
end

# === Configuration =========================================================

const TIER = Symbol(get(ENV, "BENCH_TIER", "small"))
const TRACERS = Symbol(get(ENV, "BENCH_TRACERS", "all"))

haskey(TIERS, TIER) || error("unknown BENCH_TIER=$TIER; expected one of $(collect(keys(TIERS)))")
const TRACER_KEYS = TRACERS === :all ? collect(keys(TRACER_SETUPS)) :
    haskey(TRACER_SETUPS, TRACERS) ? [TRACERS] :
    error("unknown BENCH_TRACERS=$TRACERS; expected one of $(collect(keys(TRACER_SETUPS))) or :all")

# === Termination + per-solver kwargs =======================================
# Single source of truth: `AIBECS.default_termination_condition()` returns
# `(termination_condition, abstol, reltol)` matching CTKAlg's defaults — the
# ∞-norm relative-drift criterion with `abstol = 0` and
# `reltol = 1 / (1 Myr in seconds) ≈ 3.17e-14`. Splatted into every solver
# that accepts `termination_condition`, so cross-solver comparisons share
# the same stopping rule. NLsolveJL is the lone exception: it rejects
# `termination_condition` and is L2-norm-only.
const DTC = default_termination_condition()
# 20 (not 10): CTKAlg's Shamanskii update reuses stale Jacobians and so can
# legitimately take more Newton iterations than a Jacobian-every-step
# solver while still winning on wall time, because each "free" iteration
# skips the dominant per-step cost (the sparse factorisation). Capping at
# 10 would penalise that pattern by max-itering it before convergence.
const MAXITERS = 20

const KWARGS_BY_FAMILY = Dict(
    :ctk => (; DTC..., maxItNewton = MAXITERS, preprint = "    "),
    :nr  => (; DTC..., maxiters = MAXITERS, show_trace = Val(true), trace_level = TraceMinimal()),
    # NLsolveJL: scalar abstol/reltol only, applied to the L2 norm (not
    # directly comparable to the ∞-norm criterion the other rows share).
    # The matching plot fades NLsolveJL bars to flag this.
    :nlsolve => (; abstol = DTC.abstol, reltol = DTC.reltol, maxiters = MAXITERS, show_trace = Val(true), trace_level = TraceMinimal()),
)
kwargs_for(case) = KWARGS_BY_FAMILY[case.family]

# === Stats extraction ======================================================
# Returns `(iterations, jacobian_refreshes, linear_solves)` as
# `Union{Int, Missing}`. CTKAlg + NewtonRaphson populate `sol.stats` (a
# `SciMLBase.NLStats`); NLsolveJL surfaces only `iterations` via
# `sol.original`. CTKAlg sets `nf = 0` (not tracked) — we don't read it.
function extract_stats(case, sol)
    if case.family === :ctk || case.family === :nr
        s = sol.stats
        return (iterations = s.nsteps, jacobian_refreshes = s.njacs, linear_solves = s.nsolve)
    elseif case.family === :nlsolve
        iter = sol.original === nothing ? missing : sol.original.iterations
        return (iterations = iter, jacobian_refreshes = missing, linear_solves = missing)
    else
        return (iterations = missing, jacobian_refreshes = missing, linear_solves = missing)
    end
end

# === Per-case problem builder =============================================
# Re-run between solves to clear any state the previous solve cached on the
# problem (e.g. LinearSolve symbolic factorisations, NonlinearSolve internal
# caches keyed off `prob` identity).
function fresh_problem(ssprob_unscaled, case, sp_or_nothing)
    src = sp_or_nothing === nothing ? ssprob_unscaled : sp_or_nothing.scaled
    return case.needs_nlprob ? AIBECS.nonlinearproblem(src) : src
end

# Build the `ScaledProblem` (or `nothing`) for `(setup, ssprob, scaled)`,
# reusing the unscaled `ssprob` as the source of truth.
function build_scaled(setup, ssprob; scaled::Bool)
    scaled || return nothing
    scale_u, scale_F = setup.scales(ssprob)
    return nondimensionalize(ssprob; scale_u = scale_u, scale_F = scale_F)
end

# Original-units residual + state, regardless of whether the solve ran in
# scaled or unscaled space.
function residual_orig_units(sol, sp_or_nothing, ssprob)
    if sp_or_nothing === nothing
        return sol.u, ssprob.f.f(sol.u, ssprob.p)
    else
        F̃ = sp_or_nothing.scaled.f.f(sol.u, sp_or_nothing.scaled.p)
        return sol.u .* sp_or_nothing.scale_u, F̃ .* sp_or_nothing.scale_F
    end
end

# Per-tracer implied steady-state timescale `τ = ‖u‖_v / ‖F‖_v` in *original*
# units, using a volume-weighted RMS — physically the right diagnostic for
# an ocean tracer.
function per_tracer_τ(u_orig, F_orig, v_cells, tracer_names)
    nt = length(tracer_names)
    nb_per = length(u_orig) ÷ nt
    sumv = sum(v_cells)
    vwnorm(y) = sqrt(sum(v_cells[i] * abs2(y[i]) for i in eachindex(y)) / sumv)
    return map(1:nt) do j
        idx = ((j - 1) * nb_per + 1):(j * nb_per)
        Nu = vwnorm(view(u_orig, idx))
        NF = vwnorm(view(F_orig, idx))
        NF > 0 ? Nu / NF : Inf
    end
end

# Build a failure-row NamedTuple with consistent shape — used for the skip /
# warmup-error / timed-error branches so the results array stays uniform.
function failure_row(circ, setup, case, scaled, n; retcode, error_msg = nothing)
    return (
        circulation = circ, tracers = setup.label, alg = case.label,
        scaled = scaled, n = n, seconds = NaN,
        max_rel_drift = NaN, converged = false,
        retcode = retcode,
        iterations = missing, jacobian_refreshes = missing, linear_solves = missing,
        τ_seconds = Float64[], tracer_names = setup.tracer_names,
        error = error_msg,
    )
end

# === Main loop ============================================================

vprintln("Benchmark configuration: tier=$TIER tracers=$TRACERS")
results = NamedTuple[]

for (circ_label, loader) in TIERS[TIER]
    vprintln("\n[$circ_label] loading circulation…")
    grd, T_raw = loader()
    T_circ = ustrip.(T_raw)        # toy circulations may carry units
    v_cells = volumevec(grd)
    for tk in TRACER_KEYS
        setup = TRACER_SETUPS[tk]
        vprintln("[$circ_label / $(setup.label)] building problem…")
        ssprob = setup.build(grd, T_circ)
        n = length(ssprob.u0)
        # Reference scales for the universal `max_rel_drift = max |F/scale_F|`
        # metric. When the setup ships its own physics-based `scales` we use
        # those; otherwise we fall back to the heuristic. For the unscaled
        # cells of single-tracer setups this is still a well-defined per-
        # component normalisation (≈ 1 in the surface box for idealage,
        # F0-magnitude elsewhere).
        _, ref_scale_F = setup.scales === nothing ?
            default_scales(ssprob) : setup.scales(ssprob)

        # Only run the scaled pass when the setup ships physics-based scales.
        scales_to_run = setup.scales === nothing ? (false,) : (true, false)

        for case in ALL_CASES, scaled in scales_to_run
            row_label = "$(circ_label) / $(setup.label) / $(case.label) [$(scaled ? "scaled" : "unscaled")]"
            vprintln("[$row_label] solve (n=$n)…")
            sp = build_scaled(setup, ssprob; scaled = scaled)
            kwargs = kwargs_for(case)

            # Two-pass measurement: warmup pays first-call JIT, then a fresh
            # problem is built so the timed solve sees no cached state.
            try
                solve(fresh_problem(ssprob, case, sp), case.alg; kwargs...)
            catch e
                vprintln("    ✗ solve failed (warmup): ", sprint(showerror, e))
                push!(results, failure_row(
                    circ_label, setup, case, scaled, n;
                    retcode = SciMLBase.ReturnCode.Failure,
                    error_msg = sprint(showerror, e),
                ))
                continue
            end

            timed = try
                @timed solve(fresh_problem(ssprob, case, sp), case.alg; kwargs...)
            catch e
                vprintln("    ✗ solve failed: ", sprint(showerror, e))
                push!(results, failure_row(
                    circ_label, setup, case, scaled, n;
                    retcode = SciMLBase.ReturnCode.Failure,
                    error_msg = sprint(showerror, e),
                ))
                continue
            end

            sol = timed.value
            seconds = timed.time
            u_orig, F_orig = residual_orig_units(sol, sp, ssprob)
            max_rel_drift = maximum(abs(F_orig[i]) / ref_scale_F[i] for i in eachindex(F_orig))
            converged = max_rel_drift < 1.0e-6
            stats = extract_stats(case, sol)
            τ_seconds = per_tracer_τ(u_orig, F_orig, v_cells, setup.tracer_names)
            τ_str = join(
                ("$(name)=$(@sprintf("%.2g", ustrip(t * u"s" |> u"yr"))) yr"
                    for (name, t) in zip(setup.tracer_names, τ_seconds)),
                ", ",
            )
            stats_str = "iter=$(stats.iterations)" *
                (ismissing(stats.jacobian_refreshes) ? "" : "  Jrefresh=$(stats.jacobian_refreshes)") *
                (ismissing(stats.linear_solves) ? "" : "  Lsolves=$(stats.linear_solves)")
            vprintln(
                "    → $(@sprintf("%.3g", seconds)) s",
                "   max-rel-drift = $(@sprintf("%.2e", max_rel_drift))",
                "   $stats_str   τ: $τ_str",
                "   $(converged ? "✓" : "✗")  retcode=$(sol.retcode)",
            )

            push!(results, (
                circulation = circ_label, tracers = setup.label, alg = case.label,
                scaled = scaled, n = n, seconds = seconds,
                max_rel_drift = max_rel_drift, converged = converged,
                retcode = sol.retcode,
                iterations = stats.iterations,
                jacobian_refreshes = stats.jacobian_refreshes,
                linear_solves = stats.linear_solves,
                τ_seconds = τ_seconds,
                tracer_names = setup.tracer_names,
                error = nothing,
            ))
        end
    end
end

# === Output paths =========================================================

const OUTDIR = joinpath(@__DIR__, "results")
isdir(OUTDIR) || mkpath(OUTDIR)
const MD_PATH = joinpath(OUTDIR, "latest_timings_$(TIER).md")
# Figures are written directly into the docs source tree so the next
# `docs/make.jl` run picks them up via the relative-path image references
# in `solvers.md`. One PNG per `(sweep, circ, tracer)` combination.
const FIG_DIR = normpath(joinpath(@__DIR__, "..", "docs", "src", "explanation", "figures"))
isdir(FIG_DIR) || mkpath(FIG_DIR)
fig_path(sweep, circ, tracer) = joinpath(
    FIG_DIR, "bench_$(sweep)_$(circ)_$(tracer)_$(TIER).png",
)

# === Group results ========================================================

group_order = Tuple{String, String}[]
groups = Dict{Tuple{String, String}, Vector{NamedTuple}}()
for r in results
    key = (r.circulation, r.tracers)
    if !haskey(groups, key)
        groups[key] = NamedTuple[]
        push!(group_order, key)
    end
    push!(groups[key], r)
end

fmt_int(x) = ismissing(x) ? "—" : string(x)
fmt_t(t) = isnan(t) ? "—" : (t < 1 ? @sprintf("%.2g", t) : string(round(Int, t)))
fmt_τ(t) = isinf(t) ? "∞" : @sprintf("%.2g", ustrip(t * u"s" |> u"yr"))

# Did this row's solver return a successful retcode? (Not the same as
# `converged`: a row may finish with `MaxIters` retcode but still hit a
# small `max_rel_drift`.)
row_success(r) = SciMLBase.successful_retcode(r.retcode)
rows_for_sweep(rows, sweep) = let labels = Set(c.label for c in SWEEPS[sweep].cases)
    filter(r -> r.alg in labels, rows)
end

# === Markdown table =======================================================

open(MD_PATH, "w") do io
    println(io, "<!-- Auto-generated by benchmark/solvers.jl. Do not edit by hand. -->")
    println(io)
    println(io, "# Tier: ", TIER)
    for sweep in keys(SWEEPS)
        meta = SWEEPS[sweep]
        println(io)
        println(io, "## ", meta.title, " — ", meta.subtitle)
        for (circ, tracer) in group_order
            rows = rows_for_sweep(groups[(circ, tracer)], sweep)
            isempty(rows) && continue
            n = first(rows).n
            println(io)
            println(io, "### ", circ, " / ", tracer, " (n = ", n, ")")
            println(io)
            println(io, "| Algorithm | Scale | Time (s) | Iter | Jrefresh | Lsolves | Max rel drift | τ (yr) | Retcode |")
            println(io, "|---|:---:|---:|---:|---:|---:|---:|---|:---:|")
            for r in rows
                drift_str = isnan(r.max_rel_drift) ? "—" : @sprintf("%.2e", r.max_rel_drift)
                τ_str = isempty(r.τ_seconds) ? "—" :
                    join(("$(name)=$(fmt_τ(t))"
                            for (name, t) in zip(r.tracer_names, r.τ_seconds)), ", ")
                scale_str = r.scaled ? "●" : "○"
                println(
                    io, "| ", r.alg, " | ", scale_str, " | ", fmt_t(r.seconds), " | ",
                    fmt_int(r.iterations), " | ",
                    fmt_int(r.jacobian_refreshes), " | ",
                    fmt_int(r.linear_solves), " | ",
                    drift_str, " | ", τ_str, " | ", string(r.retcode), " |",
                )
            end
        end
    end
    println(io)
    println(io, "Legend: ● = scaled (`nondimensionalize` applied); ○ = unscaled (raw original problem).")
    println(io)
    commit = get(ENV, "GITHUB_SHA", nothing)
    println(io, "_Last updated: ", string(now()),
        commit === nothing ? "" : ", commit $commit", "_")
end
vprintln("\nWrote Markdown table to ", MD_PATH)

# === Plots ================================================================
# One PNG per `(sweep, circ, tracer)`. Each figure is a 2x1 layout:
#   top    — horizontal wall-time bars (one or two per algorithm row
#            depending on `scaled` availability for this setup).
#   bottom — diagnostics scatter using the unscaled-row counts: square =
#            Newton iter, diamond = linear solves, cross = J refreshes.
# Failed-retcode rows render at alpha = 0.3 so they're visually muted.

const ALPHA_OK_SCALED   = 0.85
const ALPHA_OK_UNSCALED = 0.35
const ALPHA_FAIL        = 0.3

# Apply (ok? alpha → fail alpha) downgrade so failed rows fade out.
α(scaled, ok) = ok ? (scaled ? ALPHA_OK_SCALED : ALPHA_OK_UNSCALED) : ALPHA_FAIL

# Build one figure for the given sweep × cell.
function plot_cell(sweep, circ, tracer, rows; path)
    sweep_rows = rows_for_sweep(rows, sweep)
    isempty(sweep_rows) && return nothing
    # Algorithm ordering follows the sweep's declared order (REFAlg first).
    alg_order = [c.label for c in SWEEPS[sweep].cases]
    nalg = length(alg_order)
    n = first(rows).n

    # 2x1 stack: walltime top, diagnostics bottom. Height scales with #algs.
    fig_height = 180 + nalg * 64
    fig = Figure(size = (820, fig_height))

    ax_t = Axis(fig[1, 1];
        title = "$(SWEEPS[sweep].title) — $(circ) / $(tracer)  (n = $n)",
        xlabel = "", ylabel = "",
        yticks = (1:nalg, alg_order),
        xtickformat = xs -> [@sprintf("%.1f", x) for x in xs],
    )
    ax_d = Axis(fig[2, 1];
        xlabel = "Wall time (s)", ylabel = "",
        yticks = (1:nalg, alg_order),
        xtickformat = xs -> [@sprintf("%.1f", x) for x in xs],
    )

    # Wall-time bars.
    for (j, label) in enumerate(alg_order)
        for scaled in (true, false)
            match = filter(r -> r.alg == label && r.scaled == scaled, sweep_rows)
            isempty(match) && continue
            r = first(match)
            isnan(r.seconds) && continue
            ok = row_success(r)
            # Centre the two bars at y = j ± 0.18; if the setup only ran
            # one (unscaled), keep it centred at y = j.
            has_both = any(rr -> rr.alg == label && rr.scaled, sweep_rows) &&
                       any(rr -> rr.alg == label && !rr.scaled, sweep_rows)
            y = has_both ? j + (scaled ? -0.18 : 0.18) : j
            barplot!(
                ax_t, [r.seconds], [y];
                direction = :x, width = has_both ? 0.32 : 0.6,
                color = (:steelblue, α(scaled, ok)),
                strokecolor = :black, strokewidth = 0.5,
            )
        end
    end

    # Diagnostics scatter — unscaled rows only (or the only-row in the
    # idealage/radiocarbon case where scaled was skipped).
    for (j, label) in enumerate(alg_order)
        match = filter(r -> r.alg == label && !r.scaled, sweep_rows)
        # If unscaled isn't present (shouldn't happen with current setup)
        # fall through to the scaled row instead so the panel isn't empty.
        if isempty(match)
            match = filter(r -> r.alg == label, sweep_rows)
        end
        isempty(match) && continue
        r = first(match)
        isnan(r.seconds) && continue
        ok = row_success(r)
        a = ok ? 0.7 : ALPHA_FAIL
        for (sym, color, getter) in (
                (:rect,    :steelblue,  rr -> rr.iterations),
                (:diamond, :darkorange, rr -> rr.linear_solves),
                (:xcross,  :firebrick,  rr -> rr.jacobian_refreshes),
            )
            v = getter(r)
            ismissing(v) && continue
            scatter!(
                ax_d, [r.seconds], [j];
                marker = sym, color = (color, a), markersize = 18,
                strokewidth = 0.6, strokecolor = :black,
            )
            text!(
                ax_d, r.seconds, j;
                text = string(Int(round(v))), align = (:center, :center),
                fontsize = 9, color = (:black, a),
            )
        end
    end

    ax_t.yreversed = true
    ax_d.yreversed = true
    xlims!(ax_t, low = 0)
    xlims!(ax_d, low = 0)
    linkxaxes!(ax_t, ax_d)
    rowgap!(fig.layout, 1, 6)
    # Compact legend for both panels.
    Legend(
        fig[3, 1],
        [
            PolyElement(color = (:steelblue, ALPHA_OK_SCALED)),
            PolyElement(color = (:steelblue, ALPHA_OK_UNSCALED)),
            PolyElement(color = (:steelblue, ALPHA_FAIL)),
            MarkerElement(marker = :rect,    color = (:steelblue, 0.7), markersize = 14),
            MarkerElement(marker = :diamond, color = (:darkorange, 0.7), markersize = 14),
            MarkerElement(marker = :xcross,  color = (:firebrick, 0.7), markersize = 14),
        ],
        ["scaled", "unscaled", "failed retcode", "Newton iter", "Lin solves", "Jac refreshes"],
        orientation = :horizontal, nbanks = 2, framevisible = false,
    )
    save(path, fig)
    return path
end

if any(!isnan(r.seconds) for r in results)
    for sweep in keys(SWEEPS), (circ, tracer) in group_order
        path = fig_path(sweep, circ, tracer)
        plot_cell(sweep, circ, tracer, groups[(circ, tracer)]; path = path)
    end
    vprintln("Wrote figures to ", FIG_DIR)
else
    vprintln("No successful timings — skipping plot generation.")
end

println("\nWrote ", length(results), " benchmark entries.")
