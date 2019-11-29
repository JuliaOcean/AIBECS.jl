import AIBECS: units, @units
import AIBECS: initial_value, @initial_value
import AIBECS: flattenable, @flattenable
import AIBECS: description, @description
@flattenable @description @units @initial_value struct TestParameters{T} <: AbstractParameters{T}
    xgeo::T | 2.17 | u"mmol/m^3" | "Geological mean P concentration"             | true
    τgeo::T |  1.0 | u"Myr"      | "Geological restoring timescale"              | false
    k::T    | 10.0 | u"μmol/m^3" | "Half-saturation constant (Michaelis-Menten)" | true
    z₀::T   | 80.0 | u"m"        | "Depth of the euphotic layer base"            | false
    w₀::T   |  1.0 | u"m/d"      | "Sinking velocity at surface"                 | true
    w′::T   | 0.22 | u"m/d/m"    | "Vertical gradient of sinking velocity"       | true
    τDOP::T | 0.25 | u"yr"       | "Remineralization timescale (DOP to DIP)"     | true
    τPOP::T | 5.25 | u"d"        | "Dissolution timescale (POP to DOP)"          | true
    σ::T    |  1/3 | 1           | "Fraction of quick local uptake recycling"    | false
    τDIP::T | 30.0 | u"d"        | "Uptake maximum timescale (DIP to POP)"       | true
end

p = TestParameters(initial_value(TestParameters)...)

@testset "Parameters" begin
    println(p)
    @test p isa AIBECS.AbstractParameters{Float64}
    @testset "Adding a dual-valued vector" begin
        using DualNumbers
        dualp = p + ε * ones(sum(flattenable(p)))
        @test dualp isa AIBECS.AbstractParameters{Dual{Float64}}
        println(dualp)
    end
    @testset "Non-exported functions" begin
        @test AIBECS.symbols(p) isa NTuple{length(fieldnames(typeof(p))), Symbol}
        @test AIBECS.flattenable_symbols(p) isa NTuple{length(p), Symbol}
        @test length(values(p)) == length(fieldnames(typeof(p)))
        @test AIBECS.flattenable_values(p) isa Vector{Float64}
        @test length(AIBECS.flattenable_values(p)) == length(p)
    end
    @testset "Exported functions" begin
        @test values(p) isa Vector{Float64}
        @test vec(p) isa Vector{Float64}
        @test size(p) == size(vec(p))
        @test length(p) == length(vec(p))
        println(latex(p))
    end
    @testset "Unpacking" begin
        @unpack xgeo, τDIP = p
        @test xgeo == ustrip(upreferred(p.xgeo * units(p, :xgeo)))
        @test τDIP == ustrip(upreferred(p.τDIP * units(p, :τDIP)))
    end
end
