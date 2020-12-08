# assumes that you have run 
# - parameters.jl
# - setup.jl
# - bgc_functions.jl
# before

#
using NonlinearSolve


# C T Kelley
using SIAMFANLEquations
# p is fixed here
testF!(dx, x) = F!(dx, x, p)
nt = length(Gs) # number of tracers
x0 = xgeo * ones(nt*nb)
FS = similar(x0)
FPS = ∇ₓF(x, p)
function J!(dJ, x, p)
    dJ = ∇ₓF(x, p)
    dJ
end


nsol(F!, x0, FS, FPS, J!=diffjac!; rtol=1.e-6, atol=1.e-12,
           maxit=20, solver="newton", sham=5, armmax=10, resdec=.1,
           dx = 1.e-7, armfix=false, 
           pdata = nothing, jfact = klfact,
           printerr = true, keepsolhist = false, stagnationok=false)
