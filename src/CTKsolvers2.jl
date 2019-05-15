function updateJacobian!(JF, rShamᵢ, rSham₀, ArmijoFail, ∇ₓF, τᵢ, xᵢ; update_factors=true)
    if update_factors && rShamᵢ ≤ rSham₀ && ~ArmijoFail
        JF.age += 1
        JF.facage +=1
    else
        JF.fac = factorize(JF.fac, ∇ₓF(xᵢ) - I / τᵢ, update_factors=update_factors)
        JF.age = 0
        update_factors ? JF.facage = 0 : JF.facage += 1
    end
    return JF
end


function NewtonChordShamanskii2(F, ∇ₓF, nrm, xinit, τstop; preprint="", maxItNewton=50)
    if preprint ≠ ""
        println(preprint * "(No initial Jacobian factors fed to Newton solver)")
    end
    # Initial Jacobian
    Nx = nrm(xinit)
    NF = nrm(F(xinit))
    τ = Nx / NF
    J = ∇ₓF(xinit) - I / τ
    Jfinit = factorize(J) # construct JF
    return NewtonChordShamanskii(F, ∇ₓF, nrm, xinit, τstop, Jfinit; preprint=preprint, maxItNewton=maxItNewton)
end

function NewtonChordShamanskii2(F, ∇ₓF, nrm, xinit, τstop, Jfinit; preprint="", maxItNewton=50, update_factors=true, firstJ=false)
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

    τᵢ = Nxᵢ / NFᵢ
    rShamᵢ = 0

    # Initialize the first Jacobian
    if firstJ
        JF = AgedJacobianFactors(Jfinit, 10, 10)
    else
        JF = initialize_AgedJacobianFactor(∇ₓF(xinit) - I / τᵢ, Jfinit, update_factors)
    end

    i = 0 # counter of quasi-Newton δxⱼs
    if preprint ≠ ""
        print(preprint)
        @printf "%3s   %8s   %8s   %7s   %7s" "iteration" "|F(x)|" "|δx|/|x|" "Jac age" "fac age"
        if solType == Complex{Float64}
            @printf "%12s" "|δix|/|ix|"
        elseif solType == Dual{Float64}
            @printf "%12s" "|δεx|/|εx|"
        elseif solType == Hyper{Float64} 
            @printf "%12s" "|δhx|/|hx|"
        end
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
        JF = updateJacobian!(JF, rShamᵢ, rSham₀, ArmijoFail, ∇ₓF, τᵢ, xᵢ, update_factors=update_factors)

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
        elseif solType == Complex{Float64}
            RNimδxᵢ = nrm(imag.(δxᵢ)) / nrm(imag.(xᵢ₋₁))
            nonrealtolx = (RNimδxᵢ > 1e-7)
        elseif solType == Dual{Float64}
            RNεδxᵢ = nrm(dualpart.(δxᵢ)) / nrm(dualpart.(xᵢ₋₁))
            nonrealtolx = (RNεδxᵢ > 1e-7)
        elseif solType == Hyper{Float64}
            RNhδxᵢ = nrm(ε₁part.(δxᵢ)) / nrm(ε₁part.(xᵢ₋₁)) +
                     nrm(ε₂part.(δxᵢ)) / nrm(ε₂part.(xᵢ₋₁)) +
                     nrm(ε₁ε₂part.(δxᵢ)) / nrm(ε₁ε₂part.(xᵢ₋₁))
            nonrealtolx = (RNhδxᵢ > 1e-7)
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
            print_marker(JF.facage)
            if solType == Complex{Float64}
                @printf "%8.1e" RNimδxᵢ
            elseif solType == Dual{Float64}
                @printf "%8.1e" RNεδxᵢ
            elseif solType == Hyper{Float64} 
                @printf "%8.1e" RNhδxᵢ
            end
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

