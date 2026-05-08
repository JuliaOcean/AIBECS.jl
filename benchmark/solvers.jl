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
# Algorithm matrix: CTKAlg (baseline) plus native NonlinearSolve algorithms
# (NewtonRaphson) paired with various LinearSolve.jl factorizations. The
# wrappers SIAMFANLEquationsJL and NLsolveJL are included for cross-ecosystem
# comparison, both configured to ride the LinearSolve.jl / sparse path.
#
# NLSolversJL is intentionally omitted: its NonlinearSolve extension does not
# plumb `jac_prototype` through to NLSolvers.jl, which then allocates
# `J = similar(x, N, N)` dense (~57 GB at n≈84k). Skip until upstream lands.
#
# Krylov solvers are deliberately omitted: AIBECS lacks a preconditioner
# suitable for the steady-state advection-diffusion-reaction system.

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
import ParU_jll              # activates LinearSolveParUExt → ParUFactorization()
import SIAMFANLEquations     # activates NonlinearSolveSIAMFANLEquationsExt
import NLsolve               # activates NonlinearSolveNLsolveExt

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
# `benchmark/problems.jl`, included near the top alongside
# `dimensional_scaling.jl`. The same file is reused by `test/sparse_jacobian.jl`.

const TRACER_SETUPS = (
    # `scales = nothing` ⇒ heuristic `default_scales` from initial guess.
    # OK for single-tracer setups (no scale-mixing pathology). `po4pop`
    # uses physics-based scales because the heuristic captures `u₀`'s
    # magnitude (DIP_geo for both tracers) rather than the steady-state
    # POP magnitude (≈ DIP_geo · τ_POP/τ_DIP, ~46× smaller).
    idealage    = (label = "idealage",    build = build_idealage_problem,    scales = nothing),
    radiocarbon = (label = "radiocarbon", build = build_radiocarbon_problem, scales = nothing),
    po4pop      = (label = "po4pop",      build = build_pmodel_problem,      scales = pmodel_scales),
)

# === Algorithm matrix ======================================================

# NLsolveJL's default `linsolve = (x, A, b) -> copyto!(x, A \ b)` falls
# through to Julia's generic sparse solve. Routing the per-step solve through
# LinearSolve.jl / UMFPACK matches NewtonRaphson + UMFPACK's linear-algebra
# path, so the wrapper benchmark exercises NLsolve's nonlinear logic against
# the same factorization the native row uses.
nlsolve_linsolve_umfpack(x, A, b) =
    (copyto!(x, solve(LinearProblem(A, b), UMFPACKFactorization()).u); x)

const ALG_MATRIX = [
    # `supports_nl_kwargs = false` for wrappers (SIAMFANLEquationsJL,
    # NLsolveJL): each errors with "does not support termination conditions"
    # when `termination_condition` is passed, so they take only the universal
    # `reltol` / `abstol`.
    #
    # Both sparse-aware wrappers propagate AIBECS's analytic sparse Jacobian:
    # SIAMFANLEquationsJL via `nsol(..., FPS = zero(jac_prototype), AJ!)`, and
    # NLsolveJL via `OnceDifferentiable(f!, jac!, ..., zero(jac_prototype))`.
    # SIAMFANLEquationsJL is locked to UMFPACK by an inplace `lu!(FPS)` on the
    # prototype, which AIBECS already supplies as `SparseMatrixCSC{Float64,
    # Int64}` — the sparse path is selected automatically.
    (label = "CTKAlg",                  alg = CTKAlg(),                                         needs_nlprob = false, supports_nl_kwargs = false, max_n = Inf),
    (label = "NewtonRaphson + UMFPACK", alg = NewtonRaphson(linsolve = UMFPACKFactorization()), needs_nlprob = true,  supports_nl_kwargs = true,  max_n = Inf),
    (label = "NewtonRaphson + KLU",     alg = NewtonRaphson(linsolve = KLUFactorization()),     needs_nlprob = true,  supports_nl_kwargs = true,  max_n = Inf),
    # NewtonRaphson + ParU: known to hang inside METIS NodeND on the OCIM0
    # sparsity pattern (ParU_C_Analyze → umf_l_cholmod → cholmod_l_metis →
    # recursive MlevelNestedDissection). Kept in the matrix so we notice if /
    # when the upstream SuiteSparse / ParU_jll fix lands.
    (label = "NewtonRaphson + ParU",    alg = NewtonRaphson(linsolve = ParUFactorization()),    needs_nlprob = true,  supports_nl_kwargs = true,  max_n = Inf),
    (label = "SIAMFANLEquationsJL",     alg = SIAMFANLEquationsJL(),                            needs_nlprob = true,  supports_nl_kwargs = false, max_n = Inf),
    # NLsolveJL: :newton + LinearSolve.jl/UMFPACK adapter. :newton paired
    # with the same UMFPACK path NewtonRaphson uses isolates nonlinear-solver
    # overhead in the comparison; the default :trust_region adds J'·J /
    # Cauchy work that's not interesting here.
    (label = "NLsolveJL + Newton/UMFPACK", alg = NLsolveJL(method = :newton, linsolve = nlsolve_linsolve_umfpack), needs_nlprob = true, supports_nl_kwargs = false, max_n = Inf),
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

# CTKAlg with Inf-norm: `‖ũ‖_∞ / ‖F̃‖_∞ ≥ τstop ⇔ max-rel-drift ≤ 1/τstop`.
# `τstop = 1e6` aligns with `NL_ABSTOL = 1e-6`. `preprint = "    "` turns on
# the per-iter `i | |F(x)| | |δx|/|x| | Jac age` table.
const CTK_KWARGS = (;
    maxItNewton = 20,
    preprint = "    ",
    τstop = 1.0e6,
    nrm = Base.Fix1(maximum, abs),
)

# Per-solve diagnostic that prints, alongside the universal `max-rel-drift`,
# whatever quantity the solver's *own* internal stopping test consumes —
# evaluated on the final iterate. Reviewers see at a glance whether each
# solver's internal termination criterion agrees with the apples-to-apples
# comparison.
function format_diagnostic(case, F̃, ũ)
    Linf = maximum(abs, F̃)
    L2   = norm(F̃, 2)
    ratio = norm(ũ, Inf) / max(Linf, eps())
    return if case.label == "CTKAlg"
        @sprintf("‖ũ‖_∞ / ‖F̃‖_∞ = %.2e ≥ 1e6 ? %s",
            ratio, ratio ≥ 1e6 ? "✓" : "✗")
    elseif startswith(case.label, "NLsolveJL")
        @sprintf("‖F̃‖_2 = %.2e ≤ %.0e ? %s   [L2 (~√n × tighter than the others)]",
            L2, NL_ABSTOL, L2 ≤ NL_ABSTOL ? "✓" : "✗")
    else  # NewtonRaphson + *, SIAMFANLEquationsJL — all Inf-norm
        @sprintf("‖F̃‖_∞ = %.2e ≤ %.0e ? %s",
            Linf, NL_ABSTOL, Linf ≤ NL_ABSTOL ? "✓" : "✗")
    end
end

vprintln("Benchmark configuration: tier=$TIER tracers=$TRACERS")

results = NamedTuple[]
for (circ_label, loader) in TIERS[TIER]
    vprintln("\n[$circ_label] loading circulation…")
    grd, T_raw = loader()
    T_circ = ustrip.(T_raw)        # toy circulations may carry units
    for tk in TRACER_KEYS
        setup = TRACER_SETUPS[tk]
        vprintln("[$circ_label / $(setup.label)] building problem…")
        ssprob = setup.build(grd, T_circ)
        n = length(ssprob.u0)

        # Non-dimensionalise once per (circulation, tracer); every algorithm
        # below sees the scaled problem. `np` is the corresponding scaled
        # NonlinearProblem (with sparse `jac_prototype`) for the algorithms
        # that need one; CTKAlg consumes the SteadyStateProblem directly.
        scale_u, scale_F = setup.scales === nothing ?
            default_scales(ssprob) : setup.scales(ssprob)
        sp = nondimensionalize(ssprob; scale_u = scale_u, scale_F = scale_F)
        np = AIBECS.nonlinearproblem(sp.scaled)

        for case in ALG_MATRIX
            if n > case.max_n
                vprintln("[$circ_label / $(setup.label) / $(case.label)] skipped (n=$n > max_n=$(case.max_n); dense-Jacobian wrapper)")
                push!(
                    results, (
                        circulation = circ_label, tracers = setup.label, alg = case.label,
                        n = n, seconds = NaN,
                        max_rel_drift = NaN, converged = false,
                        error = "skipped: n>max_n",
                    )
                )
                continue
            end
            vprintln("[$circ_label / $(setup.label) / $(case.label)] solve (n=$n)…")
            alg = case.alg
            prob = case.needs_nlprob ? np : sp.scaled
            kwargs = case.supports_nl_kwargs ? NL_KWARGS_NATIVE :
                     case.needs_nlprob       ? NL_KWARGS_WRAPPER :
                                               CTK_KWARGS

            timed = try
                @timed solve(prob, alg; kwargs...)
            catch e
                vprintln("    ✗ solve failed: ", sprint(showerror, e))
                push!(
                    results, (
                        circulation = circ_label, tracers = setup.label, alg = case.label,
                        n = n, seconds = NaN,
                        max_rel_drift = NaN, converged = false,
                        error = sprint(showerror, e),
                    )
                )
                continue
            end

            sol = timed.value
            seconds = timed.time
            # Both `sol.u` and `F̃` are in scaled space (ũ ~ O(1) per
            # component, F̃ ~ relative drift per component). The universal
            # cross-solver metric `max-rel-drift = max_i |F_i / scale_F_i|`
            # is exactly `‖F̃‖_∞`.
            F̃ = sp.scaled.f.f(sol.u, sp.scaled.p)
            max_rel_drift = maximum(abs, F̃)
            converged = max_rel_drift < NL_ABSTOL
            diag = format_diagnostic(case, F̃, sol.u)
            vprintln(
                "    → $(@sprintf("%.3g", seconds)) s",
                "   max-rel-drift = $(@sprintf("%.2e", max_rel_drift))",
                "   $(converged ? "✓" : "✗")",
                "   internal: $diag",
            )

            push!(
                results, (
                    circulation = circ_label, tracers = setup.label, alg = case.label,
                    n = n,
                    seconds = seconds,
                    max_rel_drift = max_rel_drift, converged = converged,
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
    payload = [
        Dict(
            "name" => bench_name(r),
            "unit" => "s/solve",
            "value" => isnan(r.seconds) ? -1.0 : r.seconds,
            "extra" => string(
                    "n=", r.n,
                    "; max_rel_drift=", @sprintf("%.3e", r.max_rel_drift),
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
        println(io, "| Algorithm | Time (s) | Max rel drift | Converged |")
        println(io, "|---|---:|---:|:---:|")
        for r in rows
            time_str = isnan(r.seconds) ? "—" : string(round(Int, r.seconds))
            drift_str = isnan(r.max_rel_drift) ? "—" : @sprintf("%.2e", r.max_rel_drift)
            conv = r.converged ? "✓" : "✗"
            println(io, "| ", r.alg, " | ", time_str, " | ", drift_str, " | ", conv, " |")
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
