"""
    updateλs(λⱼ, λⱼ₋₁, N2Fᵢ, N2Fⱼ, N2Fⱼ₋₁)

Internal helper. Parabolic line-search update of the Armijo step size, following
Kelley (2003). Not exported.
"""
function updateλs(λⱼ, λⱼ₋₁, N2Fᵢ, N2Fⱼ, N2Fⱼ₋₁)
    # Fit a parabola for a line search of the minimum of the norm of F(xᵢ + λ δxᵢ) = p(λ)
    # modified from C.T.Kelley, 2003
    #
    # input:
    #       λⱼ = current λ
    #       λⱼ₋₁ = previous λ
    #       N2Fᵢ = |F(xᵢ)|²
    #       N2Fⱼ = |F(xᵢ + λⱼ d)|²
    #       N2Fⱼ₋₁ = |F(xᵢ + λⱼ₋₁ d)|²
    #
    # output:
    #       λⱼ₊₁, λⱼ = new values of λ's assuming a parabola
    #
    # internal parameters:
    #       σ₀ = 0.1 and σ₁ = 0.5, safeguard bounds for the linesearch

    # Set internal parameters.
    σ₀ = 0.1
    σ₁ = 0.5
    # Compute coefficients of interpolation polynomial.
    #   p(λ) = |F(xᵢ)|² + (c₁ λ + c₂ λ²) / d
    # with d = (λⱼ - λⱼ₋₁) λⱼ λⱼ₋₁ < 0.
    # If c₂ > 0 we have negative curvature and default to
    #   λⱼ₊₁ = σ₁ λ.
    # Otherwise return minimizes p(λ) but ensure that σ₀ < λⱼ₊₁ / λⱼ < σ₁
    c₂ = λⱼ₋₁ * (N2Fⱼ - N2Fᵢ) - λⱼ * (N2Fⱼ₋₁ - N2Fᵢ)
    if c₂ >= 0
        λⱼ₊₁ = σ₁ * λⱼ
        return λⱼ₊₁, λⱼ
    else
        c₁ = λⱼ^2 * (N2Fⱼ₋₁ - N2Fᵢ) - λⱼ₋₁^2 * (N2Fⱼ - N2Fᵢ)
        λⱼ₊₁ = - 0.5 * c₁ / c₂
        λⱼ₊₁ = min(max(λⱼ₊₁, σ₀ * λⱼ), σ₁ * λⱼ)
        return λⱼ₊₁, λⱼ
    end

end

"""
    searchLineArmijo!(δxᵢ, xᵢ, Fᵢ, F, maxItArmijo, nrm, preprint="")

Internal helper. In-place Armijo line search shrinking the Newton step `δxᵢ`
until ``\\|F(x + λ δx)\\|`` decreases sufficiently. Not exported.
"""
function searchLineArmijo!(δxᵢ, xᵢ, Fᵢ, F, maxItArmijo, nrm, preprint = "")
    preprint ≠ "" ? preprint = preprint * "    │" : nothing
    j = 0 # Armijo iteration counter
    α = 1.0e-4 # Armijo parameter from C. T. Kelley [2003] (K03 hereafter)
    λⱼ₋₁, λⱼ = 1.0, 1.0 # Relative step size
    # Norm of F(xᵢ) never changes
    NFᵢ = nrm(Fᵢ)
    N2Fᵢ = NFᵢ^2
    # If norm at Newton step is already good enough, or if Newton step is already to small, accept the Newton step
    if typeof(xᵢ) ≠ typeof(δxᵢ)
        xᵢ = convert(typeof(δxᵢ), xᵢ) # required for type
    end
    xⱼ = xᵢ .+ δxᵢ # Newton step
    Fⱼ = F(xⱼ)
    NFⱼ = nrm(Fⱼ)
    if (NFⱼ < (1 - α * λⱼ) * NFᵢ) || (nrm(δxᵢ) / nrm(xᵢ) < 1.0e-10) || eltype(Fᵢ) ≠ Float64
        # Update xᵢ, and Fᵢ
        xᵢ .= xⱼ
        Fᵢ .= Fⱼ
        return δxᵢ, xᵢ, Fᵢ, j, (j > maxItArmijo)
    else
        # Otherwise, start with at the middle
        λⱼ = 0.5
        xⱼ .= xᵢ .+ λⱼ .* δxᵢ
        Fⱼ .= F(xⱼ)
        NFⱼ, NFⱼ₋₁ = nrm(Fⱼ), NFⱼ
        # Initialize the squares of norms
        N2Fⱼ, N2Fⱼ₋₁ = NFⱼ^2, NFⱼ₋₁^2
        while (NFⱼ >= (1 - α * λⱼ) * NFᵢ) && (j <= maxItArmijo)
            # Update λ's using assuming a parabola (see K03)
            λⱼ, λⱼ₋₁ = updateλs(λⱼ, λⱼ₋₁, N2Fᵢ, N2Fⱼ, N2Fⱼ₋₁)
            # Update xⱼ, F(xⱼ), and norms
            xⱼ .= xᵢ .+ λⱼ .* δxᵢ
            Fⱼ .= F(xⱼ)
            NFⱼ = nrm(Fⱼ)
            if preprint ≠ ""
                j == 0 ? println("(Armijo line search)") : nothing
                print(preprint)
                RNδxⱼ = nrm(λⱼ .* δxᵢ) / nrm(xᵢ)
                @printf "%3d    %8.1e   %8.1e\n" j NFⱼ RNδxⱼ
            end
            N2Fⱼ, N2Fⱼ₋₁ = NFⱼ^2, N2Fⱼ
            j += 1
        end
        # Update δxᵢ, xᵢ, and Fᵢ
        δxᵢ .= λⱼ .* δxᵢ
        xᵢ .= xⱼ
        Fᵢ .= Fⱼ
        (preprint == "" || j ≤ maxItArmijo) ? nothing : println(preprint * "x ─> Armijo failure!")
        return δxᵢ, xᵢ, Fᵢ, j, (j > maxItArmijo)
    end
end


"""
    NewtonChordShamanskii(F, ∇ₓF, xinit, linsolve_alg, tc_cache;
                          preprint="", maxItNewton=50)
        -> (x, stats, retcode)

Newton–Chord–Shamanskii solver for `F(x) = 0`, with Armijo line search and
lazy Jacobian refresh.

Drives the solver via `F(x)` and `∇ₓF(x)`, recycling factorisations of the
Jacobian whenever the Newton residual reduces by more than `rSham₀ = 0.5`,
and refreshing otherwise. Linear solves go through a `LinearSolve.LinearCache`
built once from `linsolve_alg`; convergence is delegated to `tc_cache`, a
`NonlinearSolveBase` termination-mode cache produced by
`init(prob, mode, F(xinit), xinit; abstol, reltol)`.

This is the engine behind AIBECS's [`AIBECS.CTKAlg`](@ref) algorithm wrapper
used by `SciMLBase.solve(::SteadyStateProblem)`.

Returns the converged iterate `x`, a `stats` NamedTuple carrying
`iterations` (number of quasi-Newton steps) and `jacobian_refreshes`
(total Jacobian evaluations — the initial build plus each Shamanskii
refresh), and the `SciMLBase.ReturnCode` verdict.

Reference: C. T. Kelley (2003), *Solving Nonlinear Equations with Newton's
Method*, SIAM, Frontiers in Applied Mathematics 1.
"""
function NewtonChordShamanskii(
        F, ∇ₓF, xinit, linsolve_alg, tc_cache;
        preprint = "", maxItNewton = 50, linear_cache = nothing
    )
    if preprint ≠ ""
        println(preprint * "Solving F(x) = 0 (using Shamanskii Method)")
        preprint_end = preprint * "└─> "
        preprint = preprint * "│   "
    end

    maxItArmijo = 20 # Max iterations in Armijo line search
    rSham₀ = 0.5     # Shamanskii minimum reduction to keep old J

    # The Armijo line search needs *some* scalar magnitude; use the same norm
    # the termination mode uses if it exposes one, else default to ∞-norm.
    nrm = hasproperty(tc_cache.mode, :internalnorm) ?
        tc_cache.mode.internalnorm : Base.Fix1(maximum, abs)

    # Initial values for while loop
    Fᵢ = F(xinit)      ;  Fᵢ₋₁ = copy(Fᵢ)
    solType = eltype(Fᵢ)
    xᵢ₋₁ = convert(Vector{solType}, xinit)
    xᵢ = copy(xᵢ₋₁)
    NFᵢ = nrm(Fᵢ)   ;  NFᵢ₋₁ = NFᵢ
    δxᵢ = copy(xᵢ₋₁)

    # LinearSolve cache: built once unless the caller supplied one via
    # `linear_cache`. `linsol.A = …` flips `isfresh = true` so the next
    # `solve!` refactors; `linsol.b = …` reuses the factorization. A
    # supplied cache is used as-is — the caller is expected to have set
    # `cacheval` (the factors) and `isfresh = false` to skip the iter-1
    # factorisation; Shamanskii's ratio check + Armijo's "old J → back off
    # and refactor" branch are the safety net if the warm-start is stale.
    linsol = if linear_cache === nothing
        init(LinearProblem(∇ₓF(xinit), -copy(Fᵢ)), linsolve_alg)
    else
        linear_cache
    end
    age = 0
    # Total Jacobian evaluations done by this solve: 1 for the initial build
    # above (or the caller-supplied warm-started `linear_cache`, which we
    # still count as a starting factorisation from the solver's POV), plus
    # one per Shamanskii refresh inside the loop.
    jacobian_refreshes = 1
    rShamᵢ = 0.0
    ArmijoFail = false

    # Initial check: did xinit already converge?
    terminated = tc_cache(Fᵢ, xᵢ, xᵢ₋₁)

    i = 0 # counter of quasi-Newton δxⱼs
    if preprint ≠ ""
        print(preprint)
        @printf "%3s   %8s   %8s   %7s" "iteration" "|F(x)|" "|δx|/|x|" "Jac age"
        println("")
        print(preprint)
        @printf "%5d       %8.1e   %8s   %7s   %7s\n" i NFᵢ "" "" ""
    end

    # main iteration loop
    while !terminated && (i < maxItNewton)
        i += 1
        if preprint ≠ ""
            print(preprint)
            @printf "%5d       " i
        end
        # fill containers with current iteration values
        xᵢ₋₁ .= xᵢ
        Fᵢ₋₁ .= Fᵢ
        NFᵢ₋₁ = NFᵢ

        # Shamanskii: refresh J only if |F(xᵢ)| / |F(xᵢ₋₁)| > rSham₀ (or
        # Armijo failed last iter — see end of loop body).
        if rShamᵢ ≤ rSham₀ && !ArmijoFail
            age += 1
        else
            linsol.A = ∇ₓF(xᵢ)
            age = 0
            jacobian_refreshes += 1
        end

        # Newton step via the LinearSolve cache.
        linsol.b = -Fᵢ
        solve!(linsol)
        δxᵢ .= linsol.u

        # Update δxᵢ, xᵢ, and Fᵢ with an Armijo line search
        δxᵢ, xᵢ, Fᵢ, nArmijo, ArmijoFail = searchLineArmijo!(δxᵢ, xᵢ, Fᵢ, F, maxItArmijo, nrm, preprint)

        NFᵢ = nrm(Fᵢ)
        rShamᵢ = NFᵢ / NFᵢ₋₁

        if ArmijoFail
            if age > 0 # If Armijo fails and J is old, come back and update J
                xᵢ .= xᵢ₋₁
                Fᵢ .= Fᵢ₋₁
            else # else if J is fresh it's a complete failure
                preprint ≠ "" ? println("Complete Failure") : nothing
                return xᵢ,
                    (iterations = i, jacobian_refreshes = jacobian_refreshes),
                    SciMLBase.ReturnCode.Unstable
            end
        end

        RNδxᵢ = nrm(δxᵢ) / nrm(xᵢ₋₁)

        if preprint ≠ ""
            if ArmijoFail
                if age > 0
                    print(preprint * "Armijo Failure, but old J")
                else
                    print(preprint * "Complete Armijo failure (new J)\n")
                end
            end
            nArmijo == 0 ? nothing : print(preprint * "    └─────> ") # alignment thing
            @printf "%8.1e   %8.1e   " NFᵢ RNδxᵢ
            print_marker(age)
            println("")
        end

        # Convergence check via the NonlinearSolveBase termination cache.
        terminated = tc_cache(Fᵢ, xᵢ, xᵢ₋₁)
    end

    if preprint ≠ ""
        (terminated) && print(preprint_end * "Newton has converged, ")
        (!terminated) && print(preprint_end * "Newton has reached max iterations, ")
        @printf("‖F(x)‖ = %.2e\n", NFᵢ)
    end

    # Mirror NonlinearSolveBase: prefer the termination cache's verdict (so
    # callers see `StalledSuccess` / `Stalled` from `AbsNormSafeBest`),
    # else `MaxIters` if the Newton budget ran out.
    retcode = if terminated
        tc_cache.retcode == SciMLBase.ReturnCode.Default ?
            SciMLBase.ReturnCode.Success : tc_cache.retcode
    elseif i ≥ maxItNewton
        SciMLBase.ReturnCode.MaxIters
    else
        SciMLBase.ReturnCode.Default
    end
    stats = (iterations = i, jacobian_refreshes = jacobian_refreshes)
    return xᵢ, stats, retcode
end

function print_marker(i)
    return if i == 0
        @printf " (!)%2d    " i
    else
        @printf "%6d    " i
    end
end
