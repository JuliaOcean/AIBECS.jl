using AIBECS

wet3D, grd, T = Primeau_2x2x2.load() ;

wet3D

findall(wet3D)

iwet = findall(vec(wet3D))

nb = length(iwet)

grd

[println(box) for box in grd] ;

grd.depth_3D

grd.lat_3D

T

Matrix(T)

j = 2                            # index where we put some tracer
x = 1.0 * [i == j for i in 1:nb] # vector of 0's except for index j
x, (T * x)[j]

function G(x, p)
    τ, Ratm = p.τ, p.Ratm
    return Λ(Ratm .- x, p) - x / τ
end

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

initialize_Parameters_type(t, "C14_shoebox_parameters") # creates the type for parameters
p = C14_shoebox_parameters()                            # creates the parameters object

F, ∇ₓF = state_function_and_Jacobian(p -> T, G) # generates the state function (and its Jacobian!)
x = zeros(nb)
F(x,p)

δt = ustrip(1.0u"yr" |> u"s")
AIBECS.crank_nicolson_step(x, p, δt, F, ∇ₓF)

function time_steps(x₀, Δt, n, F, ∇ₓF)
    x_hist = [x₀]
    δt = Δt / n
    for i in 1:n
        push!(x_hist, AIBECS.crank_nicolson_step(last(x_hist), p, δt, F, ∇ₓF))
    end
    return reduce(hcat, x_hist), 0:δt:Δt
end

Δt = ustrip(7500u"yr" |> u"s") # 7500 years in seconds
x₀ = ones(5)             #
x_hist, t_hist = time_steps(x₀, Δt, 1000, F, ∇ₓF) # runs the simulation

using PyPlot
clf()
C14age_hist = -log.(x_hist) * ustrip(p.τ * u"s" |> u"yr")
plot(t_hist * ustrip(1u"s" |> u"yr"), C14age_hist')
xlabel("simulation time (years)")
ylabel("¹⁴C age (years)")
legend("box " .* string.(iwet))
title("Simulation of the evolution of ¹⁴C age with Crank-Nicolson time steps")
gcf()

prob = SteadyStateProblem(F, ∇ₓF, x, p)
x_final = solve(prob, CTKAlg()).u

C14age_final = -log.(x_final) * p.τ * u"s" .|> u"yr"
println.("box ", iwet, ": ", C14age_final);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

