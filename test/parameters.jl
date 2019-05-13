# Build the parameters type and p₀
t = empty_parameter_table()    # initialize table of parameters
add_parameter!(t, :DIPgeo, 2.3u"mmol/m^3",
    optimizable = true,
    description = "Mean DIP concentration",
    LaTeX = "\\xNUT^\\obs")
add_parameter!(t, :τg, 1.0u"Myr",
    description = "Geological restoring timescale",
    LaTeX = "\\tau_\\mathrm{geo}")
add_parameter!(t, :Umax, 1e-3u"mmol/m^3/d",
    optimizable = true,
    description = "Maximum uptake rate (Michaelis-Menten)",
    LaTeX = "U_\\mathrm{max}")
add_parameter!(t, :ku, 1e-1u"mmol/m^3",
    optimizable = true,
    description = "Half-saturation constant (Michaelis-Menten)",
    LaTeX = "k_\\vec{u}")
add_parameter!(t, :z₀, 80.0u"m",
    description = "Depth of the euphotic layer base",
    LaTeX = "z_0")
add_parameter!(t, :w₀, 1.0u"m/d",
    optimizable = true,
    description = "Sinking velocity at surface",
    LaTeX = "w_0")
add_parameter!(t, :w′, 1.0u"d^-1",
    optimizable = true,
    description = "Vertical gradient of sinking velocity",
    LaTeX = "w'")
add_parameter!(t, :κ, 0.25u"d^-1",
    optimizable = true,
    description = "Remineralization rate constant",
    LaTeX = "\\kappa")
initialize_Parameters_type(t)   # Generate the parameter type

@testset "Parameters" begin
    @test t isa DataFrames.DataFrame
    @test size(t) == (m_all, 9) # m_all is # of params
    @test size(p₀) == (m,)
    @test length(p₀) == m
    @test p₀ isa AIBECS.Parameters{Float64}
end
