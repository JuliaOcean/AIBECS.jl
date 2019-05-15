function ψtcNewtonSolve!(F, ∇ₓF, x, nrm, τstop)
  # TODO x_hist = x
    @printf "Solving F(x) = 0 (Using Pseudotransient Continuation)\n"
    nit_NS = 200
    i = 0
    s = copy(x)
    FNorm = nrm(F(x))
    FNorm₀ = copy(FNorm)
    xNorm = nrm(x)
    τ = xNorm / FNorm
    sNorm = nrm(s)
  # evaluate F at the initial iterate
  # compute the stop tolerance
  # main iteration loop
    @printf "  │ i = %d: |F(x)| = %.2g\n" i FNorm
    while (i == 0) || ((τ < τstop) && (i < nit_NS))
        i += 1
        @printf "  │ i = %d: " i
    # compute the ψtc step
        τ = xNorm / FNorm
        s .= (I / τ - ∇ₓF(x)) \ F(x)
        x .= x + s # ψtc step
        FNorm = nrm(F(x))
        xNorm = nrm(x)
        sNorm = nrm(s)
        @printf "|F(x)| = %.2g" FNorm
        @printf "; |x| = %.2g" xNorm
        @printf ", |dx|/|x| = %.2g\n" sNorm / xNorm
    end
  # on failure, set the error flag
    if (nrm(x) / nrm(F(x)) < τstop)
        @printf "  └────> Newton has reached max iterations, "
        @printf "|x|/|F(x)| = %.2g years\n" ((xNorm / FNorm) / (365 * 24 * 60 * 60))
    else    
        @printf "  └────> Newton has converged, "
        @printf "|x|/|F(x)| = %.2g years\n" ((xNorm / FNorm) / (365 * 24 * 60 * 60))
    end
  # TODO return sol, x_hist
    return x
end

export ψtcNewtonSolve!
