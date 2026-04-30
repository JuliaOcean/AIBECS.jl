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
    α = 1e-4 # Armijo parameter from C. T. Kelley [2003] (K03 hereafter)
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
    if (NFⱼ < (1 - α * λⱼ) * NFᵢ) || (nrm(δxᵢ) / nrm(xᵢ) < 1e-10) || eltype(Fᵢ) ≠ Float64
        # Update xᵢ, and Fᵢ
        xᵢ .= xⱼ
        Fᵢ .= Fⱼ
        return δxᵢ, xᵢ, Fᵢ, j, (j > maxItArmijo)
    else
        # Otherwise, start with at the middle
        λⱼ = 0.5
        xⱼ .= xᵢ .+ λⱼ .* δxᵢ
        Fⱼ .= F(xⱼ)
        NFⱼ, NFⱼ₋₁ = nrm(Fⱼ),  NFⱼ
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


function updateJacobian!(JF, rShamᵢ, rSham₀, ArmijoFail, ∇ₓF, xᵢ)
    if rShamᵢ ≤ rSham₀ && ~ArmijoFail
        JF.age += 1
    else
        JF.fac = factorize(∇ₓF(xᵢ))
        JF.age = 0
    end
    return JF
end

"""
    NewtonChordShamanskii(F, ∇ₓF, nrm, xinit, τstop; preprint="", maxItNewton=50)
    NewtonChordShamanskii(F, ∇ₓF, nrm, xinit, τstop, Jfinit; preprint="", maxItNewton=50)

Newton–Chord–Shamanskii solver for `F(x) = 0`, with Armijo line search and
lazy Jacobian refresh.

Drives the solver via `F(x)` and `∇ₓF(x)`, recycling factorisations of the
Jacobian whenever the Newton residual reduces by more than `rSham₀ = 0.5`,
and refreshing otherwise. The norm `nrm` and the stopping criterion
``\\|x\\| / \\|F(x)\\| > τ_{stop}`` follow the convention of Kelley (2003).
The second form lets the caller pass a pre-computed factorisation `Jfinit`
to skip the first Jacobian factorisation.

This is the engine behind AIBECS's [`AIBECS.CTKAlg`](@ref) algorithm wrapper
used by `SciMLBase.solve(::SteadyStateProblem)`.

Reference: C. T. Kelley (2003), *Solving Nonlinear Equations with Newton's
Method*, SIAM, Frontiers in Applied Mathematics 1.
"""
function NewtonChordShamanskii(F, ∇ₓF, nrm, xinit, τstop; preprint="", maxItNewton=50)
    if preprint ≠ ""
        println(preprint * "(No initial Jacobian factors fed to Newton solver)")
    end
    # Initial Jacobian
    J = ∇ₓF(xinit)
    Jfinit = factorize(J) # construct JF
    return NewtonChordShamanskii(F, ∇ₓF, nrm, xinit, τstop, Jfinit; preprint=preprint, maxItNewton=maxItNewton)
end

function NewtonChordShamanskii(F, ∇ₓF, nrm, xinit, τstop, Jfinit; preprint="", maxItNewton=50)
    if preprint ≠ ""
        println(preprint * "Solving F(x) = 0 (using Shamanskii Method)")
        preprint_end = preprint * "└─> "
        preprint = preprint * "│   "
    end

    maxItArmijo = 20 # Max iterations in Armijo line search
     # Max quasi Newton iterations
    rSham₀ = 0.5     # Shamanskii minimum reduction to keep old J

    # Initial values for while loop
    Fᵢ = F(xinit)      ;  Fᵢ₋₁ = copy(Fᵢ)
    solType = eltype(Fᵢ) # The type of the solution will be that of J (or f)
    xᵢ₋₁ = convert(Vector{solType}, xinit)
    xᵢ = copy(xᵢ₋₁)
    Nxᵢ = nrm(xᵢ)   ;  Nxᵢ₋₁ = Nxᵢ
    NFᵢ = nrm(Fᵢ)   ;  NFᵢ₋₁ = NFᵢ
    δxᵢ = copy(xᵢ₋₁)

    # Initialize the first Jacobian
    JF = AgedJacobianFactors(Jfinit, 0)

    τᵢ = Nxᵢ / NFᵢ
    rShamᵢ = 0

    i = 0 # counter of quasi-Newton δxⱼs
    if preprint ≠ ""
        print(preprint)
        @printf "%3s   %8s   %8s   %7s" "iteration" "|F(x)|" "|δx|/|x|" "Jac age"
        println("")
        print(preprint)
        @printf "%5d       %8.1e   %8s   %7s   %7s\n" i NFᵢ "" "" ""
    end

    # Condition on relative tolerance for imaginary part of xᵢ
    nonrealtolx = true
    ArmijoFail = false

    # main iteration loop
    while (τᵢ < τstop || nonrealtolx) && (i < maxItNewton)
        i += 1
        if preprint ≠ ""
            print(preprint)
            @printf "%5d       " i
        end
        # fill containers with current iteration values
        xᵢ₋₁ .= xᵢ
        Fᵢ₋₁ .= Fᵢ
        Nxᵢ₋₁ = Nxᵢ
        NFᵢ₋₁ = NFᵢ

        # Update Jacobian (or rather its factors, JF)
        # Shamanski: only update J if |F(xᵢ)| / |F(xᵢ₋₁)| > rSham₀
        JF = updateJacobian!(JF, rShamᵢ, rSham₀, ArmijoFail, ∇ₓF, xᵢ)

        # Newton Step
        δxᵢ .= JF.fac \ -Fᵢ

        # Update δxᵢ, xᵢ, and Fᵢ with an Armijo line search
        δxᵢ, xᵢ, Fᵢ, nArmijo, ArmijoFail = searchLineArmijo!(δxᵢ, xᵢ, Fᵢ, F, maxItArmijo, nrm, preprint)

        Nxᵢ = nrm(xᵢ)
        NFᵢ = nrm(Fᵢ)
        τᵢ = Nxᵢ / NFᵢ
        rShamᵢ = NFᵢ / NFᵢ₋₁

        if ArmijoFail
            if JF.age > 0 # If Armijo fails and J is old, come back and update J
                xᵢ .= xᵢ₋₁
                Fᵢ .= Fᵢ₋₁
            else # else if J is fresh it's a complete failure
                preprint ≠ "" ? println("Complete Failure") : nothing
                return xᵢ
            end
        end

        Nδxᵢ = nrm(δxᵢ)
        RNδxᵢ = Nδxᵢ / Nxᵢ₋₁

        # Check that the non-real part is converging
        # Here I chose to only look at the relative step size, |δx|/|x|
        if solType == Float64
            nonrealtolx = false
        end

        if preprint ≠ ""
            if ArmijoFail
                if JF.age > 0
                    print(preprint * "Armijo Failure, but old J")
                else
                    print(preprint * "Complete Armijo failure (new J): ")
                    @printf("|x|/|F(x)| = %.2g years\n", τᵢ/(365*24*60*60))
                end
            end
            nArmijo == 0 ? nothing : print(preprint * "    └─────> ") # alignment thing
            @printf "%8.1e   %8.1e   " NFᵢ RNδxᵢ
            print_marker(JF.age)
            println("")
        end
    end

    if preprint ≠ ""
        (τᵢ < τstop) && print(preprint_end * "Newton has reached max iterations, ")
        (τᵢ > τstop) && print(preprint_end * "Newton has converged, ")
        @printf("|x|/|F(x)| = %.2g years\n", τᵢ/(365*24*60*60))
    end

    return xᵢ
end

function print_marker(i)
    if i==0
        @printf " (!)%2d    " i
    else
        @printf "%6d    " i
    end
end


