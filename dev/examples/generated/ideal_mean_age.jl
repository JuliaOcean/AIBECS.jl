using AIBECS

grd, T_OCIM = AIBECS.OCIM1.load()
typeof(T_OCIM), size(T_OCIM)

T_age(p) = T_OCIM

source_age(age, p) = 1

z = vector_of_depths(grd)

minimum(z)

function sink_age(age, p)
    τ = p.τ
    return age .* (z .< 20u"m") / τ
end

sms_age(age, p) = source_age(age, p) .- sink_age(age, p)

t = empty_parameter_table()    # initialize table of parameters
add_parameter!(t, :τ, 1u"s")   # add the parameter we want (τ = 1s)
initialize_Parameters_type(t, "IdealAgeParameters")  # Generate the parameter type
t

p = IdealAgeParameters()

nb = number_of_wet_boxes(grd)  # number of wet boxes
x = ones(nb)

T_matrices = (T_age,)           # bundles all the transport matrices in a tuple
sources_minus_sinks = (sms_age,) # bundles all the source-sink functions in a tuple
F, ∇ₓF = state_function_and_Jacobian(T_matrices, sources_minus_sinks, nb) # generates the state function (and its Jacobian!)
F(x,p)

prob = SteadyStateProblem(F, ∇ₓF, x, p)

age = solve(prob, CTKAlg())

using Plots

age_in_yrs = age * u"s" .|> u"yr"

horizontalslice(age_in_yrs, grd, 2000; color=:magma)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

