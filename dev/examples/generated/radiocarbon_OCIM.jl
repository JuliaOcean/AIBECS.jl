using AIBECS

grd, T = OCIM1.load() ;

function G(x, p)
    τ, Ratm = p.τ, p.Ratm
    return Λ(Ratm .- x, p) - x / τ
end

iwet = indices_of_wet_boxes(grd)
surface_boxes = grd.depth_3D[iwet] .== grd.depth[1]

function Λ(x, p)
    λ, h = p.λ, p.h
    return λ / h * surface_boxes .* x
end

t = empty_parameter_table()                   # initialize an empty table of parameters
add_parameter!(t, :τ, 5730u"yr"/log(2)) # radioactive decay e-folding timescale
add_parameter!(t, :λ, 50u"m" / 10u"yr") # piston velocity
add_parameter!(t, :h, grd.δdepth[1])    # height of top layer
add_parameter!(t, :Ratm, 1.0u"mol/m^3") # atmospheric concentration
t

initialize_Parameters_type(t, "C14_OCIM_parameters") # creates the type for parameters
p = C14_OCIM_parameters()                            # creates the parameters object

F, ∇ₓF = state_function_and_Jacobian(p -> T, G) # generates the state function (and its Jacobian!)
nb = length(iwet)
x = zeros(nb)
F(x,p)

prob = SteadyStateProblem(F, ∇ₓF, x, p) # define the problem
R = solve(prob, CTKAlg()).u             # solve the problem

C14age = -log.(R) * p.τ * u"s" .|> u"yr"

##md # !!! warn

using Plots

horizontalslice(C14age, grd, 700; color=:viridis)

zonalaverage(C14age .|> u"hyr", grd; color=:viridis)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

