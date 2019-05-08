function ψtcNewtonSolve!(f, J, x, nrm, τstop)
  # TODO x_hist = x
    @printf "Solving f(x) = 0 (Using Pseudotransient Continuation)\n"
    nit_NS = 200
    i = 0
    s = copy(x)
    fNorm = nrm(f(x))
    fNorm₀ = copy(fNorm)
    xNorm = nrm(x)
    τ = xNorm / fNorm
    sNorm = nrm(s)
  # evaluate f at the initial iterate
  # compute the stop tolerance
  # main iteration loop
    @printf "  │ i = %d: |f(x)| = %.2g\n" i fNorm
    while (i == 0) || ((τ < τstop) && (i < nit_NS))
        i += 1
        @printf "  │ i = %d: " i
    # compute the ψtc step
        τ = xNorm / fNorm
        s .= (I / τ - J(x)) \ f(x)
        x .= x + s # ψtc step
        fNorm = nrm(f(x))
        xNorm = nrm(x)
        sNorm = nrm(s)
        @printf "|f(x)| = %.2g" fNorm
        @printf "; |x| = %.2g" xNorm
        @printf ", |dx|/|x| = %.2g\n" sNorm / xNorm
    end
  # on failure, set the error flag
    if (nrm(x) / nrm(f(x)) < τstop)
        @printf "  └────> Newton has reached max iterations, "
        @printf "|x|/|f(x)| = %.2g years\n" ((xNorm / fNorm) / (365 * 24 * 60 * 60))
    else    
        @printf "  └────> Newton has converged, "
        @printf "|x|/|f(x)| = %.2g years\n" ((xNorm / fNorm) / (365 * 24 * 60 * 60))
    end
  # TODO return sol, x_hist
    return x
end

export ψtcNewtonSolve!
