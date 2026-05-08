# Benchmark harness for AIBECS steady-state solvers.
#
# Run with:
#   julia --project=benchmark benchmark/solvers.jl
#
# Configuration via env vars:
#   BENCH_TIER    ∈ {toy, small, medium, large}            (default: small)
#   BENCH_TRACERS ∈ {idealage, radiocarbon, po4pop, all}   (default: all)
#
# Outputs in benchmark/results/:
#   - bench_results.json   : github-action-benchmark format ({name, unit, value})
#   - latest_timings.md    : human-readable Markdown table for docs inclusion
#
# Algorithm matrix: CTKAlg (baseline) plus native NonlinearSolve NewtonRaphson
# with two LinearSolve.jl factorisations (UMFPACK, KLU), and the NLsolveJL
# wrapper routed through the same LinearSolve UMFPACK path for cross-ecosystem
# comparison.
#
# Solvers intentionally omitted:
# - NLSolversJL: its NonlinearSolve extension doesn't plumb `jac_prototype`
#   through to NLSolvers.jl, which then allocates `J = similar(x, N, N)` dense
#   (~57 GB at n≈84k). Skip until upstream lands.
# - SIAMFANLEquationsJL: very slow first-call compile (~10× slower than
#   NewtonRaphson) and no per-iter trace. Removed for benchmark clarity.
# - NewtonRaphson + ParU: hangs inside METIS NodeND on OCIM0 sparsity
#   (ParU_C_Analyze → umf_l_cholmod → cholmod_l_metis → recursive
#   MlevelNestedDissection). Reproduces intermittently; removed until the
#   upstream SuiteSparse / ParU_jll fix lands.
# - Krylov: AIBECS lacks a preconditioner suitable for the steady-state
#   advection-diffusion-reaction system.

using AIBECS
using JLD2          # activates AIBECSJLD2Ext so OCIM*/OCCA loaders work
using NonlinearSolve
using NonlinearSolveBase
using LinearSolve
using SparseArrays
using JSON3
using Printf
using LinearAlgebra
using SciMLBase
using Dates
using UnPack
using Unitful
using Unitful: m, d, yr, Myr, mmol, μmol
import AIBECS: @units, units
import AIBECS: @initial_value, initial_value
import NonlinearSolveBase: AbsNormSafeBestTerminationMode
import NLsolve               # activates NonlinearSolveNLsolveExt
import Pardiso               # activates LinearSolvePardisoExt; needed for MKLPardiso* constructors

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
    small = [("OCCA", () -> OCCA.load(conservemass = true)), ("OCIM0", OCIM0.load)],
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
# `benchmark/problems.jl`, included near the top alongside
# `dimensional_scaling.jl`. The same file is reused by `test/sparse_jacobian.jl`.

const TRACER_SETUPS = (
    # `scales = nothing` ⇒ heuristic `default_scales` from initial guess.
    # OK for single-tracer setups (no scale-mixing pathology). `po4pop`
    # uses physics-based scales because the heuristic captures `u₀`'s
    # magnitude (DIP_geo for both tracers) rather than the steady-state
    # POP magnitude (≈ DIP_geo · τ_POP/τ_DIP, ~46× smaller).
    # `tracer_names` is used to label per-tracer τ = ‖x‖_∞ / ‖F‖_∞ in the
    # output (see post-solve block below). The vector's length determines how
    # the state vector is sliced into per-tracer blocks (length(u) ÷ nt).
    idealage    = (label = "idealage",    build = build_idealage_problem,    scales = nothing,       tracer_names = ["age"]),
    radiocarbon = (label = "radiocarbon", build = build_radiocarbon_problem, scales = nothing,       tracer_names = ["R"]),
    po4pop      = (label = "po4pop",      build = build_pmodel_problem,      scales = pmodel_scales, tracer_names = ["DIP", "POP"]),
)

# === Algorithm matrix ======================================================

# NLsolveJL's default `linsolve = (x, A, b) -> copyto!(x, A \ b)` falls
# through to Julia's generic sparse solve. Routing the per-step solve through
# LinearSolve.jl / UMFPACK matches NewtonRaphson + UMFPACK's linear-algebra
# path, so the wrapper benchmark exercises NLsolve's nonlinear logic against
# the same factorization the native row uses.
nlsolve_linsolve_umfpack(x, A, b) =
    (copyto!(x, solve(LinearProblem(A, b), UMFPACKFactorization()).u); x)

# MKL Pardiso rows are platform-conditional. `Pardiso.mkl_is_available()`
# returns true on Linux (where MKL_jll ships native binaries) and false on
# macOS, so dev runs on a Mac silently skip the 4 Pardiso rows. `nprocs` is
# the OMP-thread count for Pardiso; on `ubuntu-latest` runners
# `Sys.CPU_THREADS` picks up the runner's full core count. `BENCH_NPROCS`
# overrides if a different value is wanted.
const PARDISO_AVAILABLE = Pardiso.mkl_is_available()
const PARDISO_NPROCS = parse(Int, get(ENV, "BENCH_NPROCS", string(Sys.CPU_THREADS)))
if PARDISO_AVAILABLE
    vprintln("MKL Pardiso available — adding 4 rows with nprocs = $PARDISO_NPROCS")
else
    vprintln("MKL Pardiso not available on this platform — skipping 4 rows")
end

const ALG_MATRIX = [
    # `supports_nl_kwargs = false` for the NLsolveJL wrapper: it errors with
    # "does not support termination conditions" when `termination_condition`
    # is passed, so it takes only the universal `reltol` / `abstol`. The
    # wrapper propagates AIBECS's analytic sparse Jacobian via
    # `OnceDifferentiable(f!, jac!, ..., zero(jac_prototype))`.
    #
    # NLsolveJL: `:newton` + LinearSolve.jl/UMFPACK adapter. The UMFPACK path
    # matches NewtonRaphson + UMFPACK so the row isolates nonlinear-solver
    # overhead in the comparison; the default `:trust_region` adds J'·J /
    # Cauchy work that's not interesting here.
    (label = "CTKAlg + UMFPACK",        alg = CTKAlg(linsolve = UMFPACKFactorization()),        needs_nlprob = false, supports_nl_kwargs = false, max_n = Inf),
    (label = "CTKAlg + KLU",            alg = CTKAlg(linsolve = KLUFactorization()),            needs_nlprob = false, supports_nl_kwargs = false, max_n = Inf),
    (label = "NewtonRaphson + UMFPACK", alg = NewtonRaphson(linsolve = UMFPACKFactorization()), needs_nlprob = true,  supports_nl_kwargs = true,  max_n = Inf),
    (label = "NewtonRaphson + KLU",     alg = NewtonRaphson(linsolve = KLUFactorization()),     needs_nlprob = true,  supports_nl_kwargs = true,  max_n = Inf),
    (label = "NLsolveJL + Newton/UMFPACK", alg = NLsolveJL(method = :newton, linsolve = nlsolve_linsolve_umfpack), needs_nlprob = true, supports_nl_kwargs = false, max_n = Inf),
    # MKL Pardiso × {NewtonRaphson, CTKAlg} — appended only when MKL Pardiso
    # is installed. Both `MKLPardisoFactorize` and `MKLPardisoIterate` accept
    # an `nprocs` kwarg (LinearSolve src/extension_algs.jl:137,165) that
    # configures Pardiso's OMP parallelism.
    (PARDISO_AVAILABLE ? [
        (label = "NewtonRaphson + MKLPardisoFactorize", alg = NewtonRaphson(linsolve = MKLPardisoFactorize(; nprocs = PARDISO_NPROCS)), needs_nlprob = true,  supports_nl_kwargs = true,  max_n = Inf),
        (label = "NewtonRaphson + MKLPardisoIterate",   alg = NewtonRaphson(linsolve = MKLPardisoIterate(; nprocs = PARDISO_NPROCS)),   needs_nlprob = true,  supports_nl_kwargs = true,  max_n = Inf),
        (label = "CTKAlg + MKLPardisoFactorize",        alg = CTKAlg(linsolve = MKLPardisoFactorize(; nprocs = PARDISO_NPROCS)),        needs_nlprob = false, supports_nl_kwargs = false, max_n = Inf),
        (label = "CTKAlg + MKLPardisoIterate",          alg = CTKAlg(linsolve = MKLPardisoIterate(; nprocs = PARDISO_NPROCS)),          needs_nlprob = false, supports_nl_kwargs = false, max_n = Inf),
    ] : [])...,
]

# === Configuration =========================================================

const TIER = Symbol(get(ENV, "BENCH_TIER", "small"))
const TRACERS = Symbol(get(ENV, "BENCH_TRACERS", "all"))

haskey(TIERS, TIER) || error("unknown BENCH_TIER=$TIER; expected one of $(collect(keys(TIERS)))")
const TRACER_KEYS = TRACERS === :all ? collect(keys(TRACER_SETUPS)) :
    haskey(TRACER_SETUPS, TRACERS) ? [TRACERS] :
    error("unknown BENCH_TRACERS=$TRACERS; expected one of $(collect(keys(TRACER_SETUPS))) or :all")

# === Run benchmarks ========================================================

# After non-dim, both `ũ` and `F̃` are O(1) per component, so a scalar
# `abstol` on the Inf-norm of `F̃` reads literally as "max relative drift per
# component ≤ NL_ABSTOL". Each algorithm's internal norm is aligned to Inf
# where configurable; NLsolveJL is stuck at L2 because the SciML extension
# (`NonlinearSolveNLsolveExt.jl:58-62`) doesn't pass a `norm` kwarg through to
# `NLsolve.nlsolve`. Apples-to-apples comparison is done post-solve via
# `max_i |F_i / scale_F_i|`, which is reported in the MD `Max rel drift`
# column.
const NL_ABSTOL = 1.0e-6

const NL_KWARGS_NATIVE = (;
    abstol = NL_ABSTOL,
    maxiters = 20,
    show_trace = Val(true),
    trace_level = TraceMinimal(),
    termination_condition = AbsNormSafeBestTerminationMode(
        Base.Fix1(maximum, abs); max_stalled_steps = 8),
)

const NL_KWARGS_WRAPPER = (;
    abstol = NL_ABSTOL,
    maxiters = 20,
    show_trace = Val(true),
    trace_level = TraceMinimal(),
)

# CTKAlg now goes through NonlinearSolveBase termination, so the abstol +
# termination_condition match `NL_KWARGS_NATIVE` exactly — the CTK row is
# truly apples-to-apples with NewtonRaphson on the residual ∞-norm.
# CTKAlg's iteration cap is `maxItNewton` (not `maxiters`) and its trace
# knob is `preprint` (not `show_trace`/`trace_level`); otherwise identical.
const CTK_KWARGS = (;
    abstol = NL_ABSTOL,
    maxItNewton = 20,
    preprint = "    ",
    termination_condition = AbsNormSafeBestTerminationMode(
        Base.Fix1(maximum, abs); max_stalled_steps = 8),
)

# Per-solve diagnostic that prints, alongside the universal `max-rel-drift`,
# whatever quantity the solver's *own* internal stopping test consumes —
# evaluated on the final iterate. Reviewers see at a glance whether each
# solver's internal termination criterion agrees with the apples-to-apples
# comparison. CTKAlg now uses the same `AbsNormSafeBestTerminationMode` as
# the NewtonRaphson rows, so its diagnostic collapses into the same branch.
function format_diagnostic(case, F̃)
    Linf = maximum(abs, F̃)
    L2   = norm(F̃, 2)
    return if startswith(case.label, "NLsolveJL")
        @sprintf("‖F̃‖_2 = %.2e ≤ %.0e ? %s   [L2 (~√n × tighter than the others)]",
            L2, NL_ABSTOL, L2 ≤ NL_ABSTOL ? "✓" : "✗")
    else  # NewtonRaphson + * and CTKAlg + * — Inf-norm + AbsNormSafeBest
        @sprintf("‖F̃‖_∞ = %.2e ≤ %.0e ? %s",
            Linf, NL_ABSTOL, Linf ≤ NL_ABSTOL ? "✓" : "✗")
    end
end

vprintln("Benchmark configuration: tier=$TIER tracers=$TRACERS")

results = NamedTuple[]

# Build a fresh `(ScaledProblem, prob_for_solve)` pair for the given case.
# Re-run between solves to clear any state the previous solve cached on the
# problem (e.g. LinearSolve symbolic factorisations, NonlinearSolve internal
# caches keyed off `prob` identity). CTKAlg consumes the SteadyStateProblem
# directly; everyone else gets the scaled NonlinearProblem.
function fresh_pair(setup, case, grd, T_circ)
    ssprob_local = setup.build(grd, T_circ)
    scale_u_local, scale_F_local = setup.scales === nothing ?
        default_scales(ssprob_local) : setup.scales(ssprob_local)
    sp_local = nondimensionalize(
        ssprob_local; scale_u = scale_u_local, scale_F = scale_F_local
    )
    prob_local = case.needs_nlprob ? AIBECS.nonlinearproblem(sp_local.scaled) :
        sp_local.scaled
    return sp_local, prob_local
end

for (circ_label, loader) in TIERS[TIER]
    vprintln("\n[$circ_label] loading circulation…")
    grd, T_raw = loader()
    T_circ = ustrip.(T_raw)        # toy circulations may carry units
    # Cell volumes for the volume-weighted 2-norm used in the per-tracer τ.
    v_cells = volumevec(grd)
    for tk in TRACER_KEYS
        setup = TRACER_SETUPS[tk]
        vprintln("[$circ_label / $(setup.label)] building problem…")
        # `n` is needed for the per-case skip check and the MD section heading;
        # the actual problems used by `solve` are rebuilt fresh per algorithm.
        n = length(setup.build(grd, T_circ).u0)

        for case in ALG_MATRIX
            if n > case.max_n
                vprintln("[$circ_label / $(setup.label) / $(case.label)] skipped (n=$n > max_n=$(case.max_n); dense-Jacobian wrapper)")
                push!(
                    results, (
                        circulation = circ_label, tracers = setup.label, alg = case.label,
                        n = n, seconds = NaN,
                        max_rel_drift = NaN, converged = false,
                        τ_seconds = Float64[], tracer_names = setup.tracer_names,
                        error = "skipped: n>max_n",
                    )
                )
                continue
            end
            vprintln("[$circ_label / $(setup.label) / $(case.label)] solve (n=$n)…")
            alg = case.alg
            kwargs = case.supports_nl_kwargs ? NL_KWARGS_NATIVE :
                     case.needs_nlprob       ? NL_KWARGS_WRAPPER :
                                               CTK_KWARGS

            # Two-pass measurement per `(circ, tracer, alg)` row:
            #   1. Warmup solve on a fresh problem — pays first-call JIT
            #      compilation so the timed solve measures actual solve cost
            #      rather than `solve + compile`. Untimed.
            #   2. Build a new fresh problem (so any state the warmup left on
            #      the previous `prob` is cleared) and `@timed solve` it.
            sp_warm, prob_warm = fresh_pair(setup, case, grd, T_circ)
            try
                solve(prob_warm, alg; kwargs...)
            catch e
                vprintln("    ✗ solve failed (warmup): ", sprint(showerror, e))
                push!(
                    results, (
                        circulation = circ_label, tracers = setup.label, alg = case.label,
                        n = n, seconds = NaN,
                        max_rel_drift = NaN, converged = false,
                        τ_seconds = Float64[], tracer_names = setup.tracer_names,
                        error = sprint(showerror, e),
                    )
                )
                continue
            end

            sp_fresh, prob_fresh = fresh_pair(setup, case, grd, T_circ)
            timed = try
                @timed solve(prob_fresh, alg; kwargs...)
            catch e
                vprintln("    ✗ solve failed: ", sprint(showerror, e))
                push!(
                    results, (
                        circulation = circ_label, tracers = setup.label, alg = case.label,
                        n = n, seconds = NaN,
                        max_rel_drift = NaN, converged = false,
                        τ_seconds = Float64[], tracer_names = setup.tracer_names,
                        error = sprint(showerror, e),
                    )
                )
                continue
            end

            sol = timed.value
            seconds = timed.time
            # `sol.u` and `F̃` are in scaled space (ũ ~ O(1) per component, F̃ ~
            # relative drift per component). The universal cross-solver metric
            # `max-rel-drift = max_i |F_i / scale_F_i|` is exactly `‖F̃‖_∞`.
            F̃ = sp_fresh.scaled.f.f(sol.u, sp_fresh.scaled.p)
            max_rel_drift = maximum(abs, F̃)
            converged = max_rel_drift < NL_ABSTOL
            diag = format_diagnostic(case, F̃)
            # Per-tracer implied steady-state timescale, computed in *original*
            # (unscaled) units so the result is genuinely in seconds (since F
            # is du/dt) and converts to years cleanly. For each tracer block,
            # `τ = ‖x‖_v / ‖F‖_v` where `‖y‖_v = √(yᵀ Diagonal(v) y / Σv)` is
            # the volume-weighted RMS — physically the right diagnostic for
            # an ocean tracer (cells with more water count more). Used only
            # for *display*; the convergence test still uses the unweighted
            # Inf-norm of the scaled residual. Values "Inf" mean the tracer's
            # residual is at floating-point zero.
            u_orig = sol.u .* sp_fresh.scale_u
            F_orig = F̃ .* sp_fresh.scale_F
            nt = length(setup.tracer_names)
            nb_per = length(u_orig) ÷ nt
            sumv = sum(v_cells)
            vwnorm(y, v) = sqrt(sum(v[i] * abs2(y[i]) for i in eachindex(y)) / sumv)
            τ_seconds = map(1:nt) do j
                idx = ((j - 1) * nb_per + 1):(j * nb_per)
                Nu = vwnorm(view(u_orig, idx), v_cells)
                NF = vwnorm(view(F_orig, idx), v_cells)
                NF > 0 ? Nu / NF : Inf
            end
            τ_str = join(
                ("$(name) = $(@sprintf("%.2g", ustrip(t * u"s" |> u"yr"))) yr"
                    for (name, t) in zip(setup.tracer_names, τ_seconds)),
                ", ",
            )
            vprintln(
                "    → $(@sprintf("%.3g", seconds)) s",
                "   max-rel-drift = $(@sprintf("%.2e", max_rel_drift))",
                "   τ: $τ_str",
                "   $(converged ? "✓" : "✗")",
                "   internal: $diag",
            )

            push!(
                results, (
                    circulation = circ_label, tracers = setup.label, alg = case.label,
                    n = n,
                    seconds = seconds,
                    max_rel_drift = max_rel_drift, converged = converged,
                    τ_seconds = τ_seconds,
                    tracer_names = setup.tracer_names,
                    error = nothing,
                )
            )
        end
    end
end

# === Write outputs =========================================================

bench_name(r) = string(r.circulation, " / ", r.tracers, " / ", r.alg)

const OUTDIR = joinpath(@__DIR__, "results")
isdir(OUTDIR) || mkpath(OUTDIR)

# One pair of files per tier so successive runs of different tiers
# accumulate on disk instead of clobbering each other.
const JSON_PATH = joinpath(OUTDIR, "bench_results_$(TIER).json")
const MD_PATH = joinpath(OUTDIR, "latest_timings_$(TIER).md")

open(JSON_PATH, "w") do io
    # `value` is what `github-action-benchmark` plots as the dashboard time
    # series. Each measurement is a single timed solve preceded by an
    # untimed warmup (so JIT compile is amortised).
    payload = [
        Dict(
            "name" => bench_name(r),
            "unit" => "s/solve",
            "value" => isnan(r.seconds) ? -1.0 : r.seconds,
            "extra" => string(
                    "n=", r.n,
                    "; max_rel_drift=", @sprintf("%.3e", r.max_rel_drift),
                    isempty(r.τ_seconds) ? "" : string(
                        "; τ_yr=",
                        join(
                            ("$(name)=$(@sprintf("%.3g", ustrip(t * u"s" |> u"yr")))"
                                for (name, t) in zip(r.tracer_names, r.τ_seconds)),
                            ",",
                        ),
                    ),
                    r.converged ? "" : "; CONVERGENCE FAILED",
                    r.error === nothing ? "" : string("; error=", r.error),
                ),
        ) for r in results
    ]
    JSON3.pretty(io, payload)
end

# One table per (circulation, tracers) cell of the run matrix — a head-to-head
# algorithm comparison where `n` is constant within the group, so it lifts
# into the section heading rather than a column.
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

open(MD_PATH, "w") do io
    println(io, "<!-- Auto-generated by benchmark/solvers.jl. Do not edit by hand. -->")
    println(io)
    println(io, "# Tier: ", TIER)
    for (circ, tracer) in group_order
        rows = groups[(circ, tracer)]
        n = first(rows).n
        println(io)
        println(io, "## ", circ, " / ", tracer, " (n = ", n, ")")
        println(io)
        # `Time (s)` is a single timed solve, preceded by an untimed warmup
        # so JIT compile cost isn't included. Sub-second solves use 2-sig-fig
        # `%.2g`; ≥1 s rounds to nearest second. `τ (yr)` is the per-tracer
        # implied steady-state timescale `‖x‖_v / ‖F‖_v` (volume-weighted RMS
        # in *original* units, so the year unit is meaningful).
        println(io, "| Algorithm | Time (s) | Max rel drift | τ (yr) | Converged |")
        println(io, "|---|---:|---:|---|:---:|")
        fmt_t(t) = isnan(t) ? "—" : (t < 1 ? @sprintf("%.2g", t) : string(round(Int, t)))
        fmt_τ(t) = isinf(t) ? "∞" : @sprintf("%.2g", ustrip(t * u"s" |> u"yr"))
        for r in rows
            time_str = fmt_t(r.seconds)
            drift_str = isnan(r.max_rel_drift) ? "—" : @sprintf("%.2e", r.max_rel_drift)
            τ_str = isempty(r.τ_seconds) ? "—" :
                join(
                    ("$(name)=$(fmt_τ(t))"
                        for (name, t) in zip(r.tracer_names, r.τ_seconds)),
                    ", ",
                )
            conv = r.converged ? "✓" : "✗"
            println(io, "| ", r.alg, " | ", time_str, " | ", drift_str, " | ", τ_str, " | ", conv, " |")
        end
    end
    println(io)
    commit = get(ENV, "GITHUB_SHA", nothing)
    println(
        io, "_Last updated: ", string(now()),
        commit === nothing ? "" : ", commit $commit", "_"
    )
end

println("\nWrote ", length(results), " benchmark entries to ", OUTDIR)
