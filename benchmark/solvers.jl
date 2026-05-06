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
import NonlinearSolveBase: RelNormSafeBestTerminationMode
import ParU_jll              # activates LinearSolveParUExt → ParUFactorization()
import SIAMFANLEquations     # activates NonlinearSolveSIAMFANLEquationsExt
import NLsolve               # activates NonlinearSolveNLsolveExt

# Match CTKAlg's default τstop (the threshold on |x| / |F(x)| in seconds).
# With reltol = 1/τstop and the Rel*Norm termination mode, NonlinearSolve's
# convergence test reduces near the root to |F(u)| ≤ |u| / τstop, which is
# the same scale-aware criterion CTKAlg applies. CTKAlg uses an L2 norm,
# so we pass `LinearAlgebra.norm` as the internal norm.
const τstop_seconds = ustrip(u"s", 1.0e6u"Myr")

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

# Ideal age — mirrors docs/lit/tutorials/1_ideal_age.jl
struct IdealAgeParameters{U} <: AbstractParameters{U}
    τ::U
    z₀::U
end
function build_idealage_problem(grd, T_circ)
    z = depthvec(grd)
    G(x, p) = (@unpack τ, z₀ = p; @. 1 - x / τ * (z ≤ z₀))
    F = AIBECSFunction(T_circ, G)
    p = IdealAgeParameters(1.0, z[1])
    nb = sum(iswet(grd))
    return SteadyStateProblem(F, zeros(nb), p)
end

# Radiocarbon — mirrors docs/lit/tutorials/2_radiocarbon.jl
@units struct RadiocarbonParameters{U} <: AbstractParameters{U}
    λ::U | u"m/yr"
    h::U | u"m"
    τ::U | u"yr"
    Ratm::U | u"M"
end
function build_radiocarbon_problem(grd, T_circ)
    z = depthvec(grd)
    G(R, p) = (@unpack λ, h, Ratm, τ = p; @. λ / h * (Ratm - R) * (z ≤ h) - R / τ)
    F = AIBECSFunction(T_circ, G)
    p = RadiocarbonParameters(
        λ = 50u"m" / 10u"yr",
        h = grd.δdepth[1],
        τ = 5730u"yr" / log(2),
        Ratm = 42.0u"nM",
    )
    nb = sum(iswet(grd))
    return SteadyStateProblem(F, zeros(nb), p)
end

# Coupled PO4–POP — mirrors docs/lit/tutorials/3_Pmodel.jl
@initial_value @units struct PmodelParameters{U} <: AbstractParameters{U}
    w₀::U | 0.64 | m / d
    w′::U | 0.13 | m / d / m
    τ_DIP::U | 230.0 | d
    k::U | 6.62 | μmol / m^3
    z₀::U | 80.0 | m
    τ_POP::U | 5.0 | d
    τ_geo::U | 1.0 | Myr
    DIP_geo::U | 2.12 | mmol / m^3
end
function build_pmodel_problem(grd, T_circ)
    nb = sum(iswet(grd))
    z = depthvec(grd)
    T_DIP(p) = T_circ
    function w(z, p)
        @unpack w₀, w′ = p
        return @. w₀ + w′ * z
    end
    T_POP(p) = transportoperator(grd, z -> w(z, p))
    function U(x, p)
        @unpack τ_DIP, k, z₀ = p
        return @. x / τ_DIP * x / (x + k) * (z ≤ z₀) * (x ≥ 0)
    end
    function R(x, p)
        @unpack τ_POP = p
        return x / τ_POP
    end
    function G_DIP(DIP, POP, p)
        @unpack DIP_geo, τ_geo = p
        return @. -$U(DIP, p) + $R(POP, p) + (DIP_geo - DIP) / τ_geo
    end
    function G_POP(DIP, POP, p)
        @unpack τ_geo = p
        return @. $U(DIP, p) - $R(POP, p) - POP / τ_geo
    end
    F = AIBECSFunction((T_DIP, T_POP), (G_DIP, G_POP), nb)
    p = PmodelParameters()
    @unpack DIP_geo = p
    return SteadyStateProblem(F, DIP_geo * ones(2nb), p)
end

const TRACER_SETUPS = (
    idealage = (label = "idealage", build = build_idealage_problem),
    radiocarbon = (label = "radiocarbon", build = build_radiocarbon_problem),
    po4pop = (label = "po4pop", build = build_pmodel_problem),
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

# Native NonlinearSolve algorithms: the Rel*Norm termination_condition with
# reltol = 1/τstop reproduces CTKAlg's scale-aware |F(u)| ≤ |u| / τstop test.
# `abstol = eps(Float64)` is harmless here — RelNormSafeBest terminates by
# stagnation (`max_stalled_steps`) once the relative criterion plateaus.
# CTKAlg itself accepts no kwargs from this surface (its API is τstop / nrm /
# maxItNewton).
const NL_TOLS = (;
    reltol = 1 / τstop_seconds,
    abstol = eps(Float64),
)
const NL_KWARGS = (;
    NL_TOLS...,
    show_trace = Val(true),
    trace_level = TraceMinimal(),
    termination_condition = RelNormSafeBestTerminationMode(
        LinearAlgebra.norm; max_stalled_steps = 8),
)

# Wrapper algorithms (SIAMFANLEquationsJL, NLsolveJL) reject
# `termination_condition` and treat `abstol` as the operative |F| tolerance.
# `eps(Float64)` is far below the residual floor at n≈10⁵ (matvec roundoff
# alone is orders of magnitude larger), so leaving it in would push NLsolve's
# `:newton` to its 1000-iter cap, silently re-factorising the same constant
# Jacobian on every step. Use 1e-8 (NLsolve's own default ftol) and cap
# iterations defensively. `show_trace` keeps these wrappers as visible as
# the native rows.
const NL_WRAPPER_KWARGS = (;
    abstol = 1.0e-8,
    maxiters = 50,
    show_trace = Val(true),
    trace_level = TraceMinimal(),
)

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
        for case in ALG_MATRIX
            if n > case.max_n
                vprintln("[$circ_label / $(setup.label) / $(case.label)] skipped (n=$n > max_n=$(case.max_n); dense-Jacobian wrapper)")
                push!(
                    results, (
                        circulation = circ_label, tracers = setup.label, alg = case.label,
                        n = n, seconds = NaN,
                        residual = NaN, converged = false,
                        error = "skipped: n>max_n",
                    )
                )
                continue
            end
            vprintln("[$circ_label / $(setup.label) / $(case.label)] solve (n=$n)…")
            alg = case.alg
            prob = case.needs_nlprob ? AIBECS.nonlinearproblem(ssprob) : ssprob
            kwargs = case.supports_nl_kwargs ? NL_KWARGS :
                     case.needs_nlprob ? NL_WRAPPER_KWARGS : (;)

            timed = try
                @timed solve(prob, alg; kwargs...)
            catch e
                vprintln("    ✗ solve failed: ", sprint(showerror, e))
                push!(
                    results, (
                        circulation = circ_label, tracers = setup.label, alg = case.label,
                        n = n, seconds = NaN,
                        residual = NaN, converged = false,
                        error = sprint(showerror, e),
                    )
                )
                continue
            end

            sol = timed.value
            seconds = timed.time
            residual = norm(ssprob.f.f(sol.u, ssprob.p))
            converged = residual < 1.0e-6 * max(norm(sol.u), 1.0)
            vprintln("    → $(@sprintf("%.3g", seconds)) s   residual=$(@sprintf("%.2e", residual))   $(converged ? "✓" : "✗")")

            push!(
                results, (
                    circulation = circ_label, tracers = setup.label, alg = case.label,
                    n = n,
                    seconds = seconds,
                    residual = residual, converged = converged,
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
                    "; residual=", @sprintf("%.3e", r.residual),
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
        println(io, "| Algorithm | Time (s) | Residual | Converged |")
        println(io, "|---|---:|---|:---:|")
        for r in rows
            time_str = isnan(r.seconds) ? "—" : string(round(Int, r.seconds))
            res_str = isnan(r.residual) ? "—" : @sprintf("%.2e", r.residual)
            conv = r.converged ? "✓" : "✗"
            println(io, "| ", r.alg, " | ", time_str, " | ", res_str, " | ", conv, " |")
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
