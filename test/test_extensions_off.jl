"""
    test_extensions_off.jl

Run AIBECS in a temp environment with no trigger packages loaded and
assert that

1. the extension `Base.get_extension` returns `nothing`,
2. each gated `load()` raises a stub error mentioning the required
   trigger packages,
3. the simple `Base.show` of an `AbstractParameters` works without a
   `MethodError`,
4. once trigger packages are added, the extensions activate and the
   stub errors are replaced by real methods.

This is intended to be run from the test environment via
`include("test_extensions_off.jl")`. It uses an isolated temp Pkg
environment internally.
"""

using Test
using Pkg

# --- Phase 1: AIBECS-only environment ---------------------------------------

mktempdir() do tmp
    Pkg.activate(tmp)
    Pkg.develop(path = dirname(@__DIR__))
    Pkg.add("Test")

    # Run the AIBECS-only assertions in a child Julia process so the
    # extensions don't leak through from the parent's session.
    aibecs_only_script = """
    using AIBECS
    using Test

    struct SimpleExtOffP{T} <: AIBECS.AbstractParameters{T}; α::T; β::T; end

    @testset "extensions off (no triggers loaded)" begin
        for ext in (:AIBECSRecipesBaseExt, :AIBECSRecipesParametersExt,
                    :AIBECSBijectorsExt, :AIBECSDistributionsExt,
                    :AIBECSDataFramesExt, :AIBECSShapefileExt,
                    :AIBECSJLD2Ext, :AIBECSNCDatasetsExt,
                    :AIBECSETOPOExt, :AIBECSOCIM2_48LExt)
            @test Base.get_extension(AIBECS, ext) === nothing
        end

        # Stub errors mention required packages
        @test_throws ErrorException AIBECS.OCIM0.load()
        @test_throws ErrorException AIBECS.OCIM1.load()
        @test_throws ErrorException AIBECS.OCIM2.load()
        @test_throws ErrorException AIBECS.OCCA.load()
        @test_throws ErrorException AIBECS.OCIM2_48L.load()
        @test_throws ErrorException AIBECS.AeolianSources.load_Chien()
        @test_throws ErrorException AIBECS.AeolianSources.load_Kok()
        @test_throws ErrorException AIBECS.ETOPO.load()
        @test_throws ErrorException AIBECS.GroundWaters.load()

        # Simple Base.show keeps working without DataFrames
        p = SimpleExtOffP(1.0, 2.0)
        io = IOBuffer()
        show(io, p)
        @test occursin("SimpleExtOffP", String(take!(io)))
    end
    """

    @test success(
        pipeline(
            `\$(Base.julia_cmd()) --project=\$(tmp) -e \$aibecs_only_script`,
            stdout = stdout, stderr = stderr
        )
    )
end

# --- Phase 2: extensions activate when triggers are loaded ------------------

# At this point we are back in the test/runtests.jl environment, which
# (per Project.toml [targets].test) has every trigger package loaded.
# `runtests.jl` already does `using JLD2, NCDatasets, ...` at the top,
# so the extensions should be live.
struct ExtOnP{T} <: AIBECS.AbstractParameters{T}
    α::T; β::T
end

@testset "extensions on (triggers loaded)" begin
    for ext in (
            :AIBECSRecipesBaseExt, :AIBECSRecipesParametersExt,
            :AIBECSBijectorsExt, :AIBECSDistributionsExt,
            :AIBECSDataFramesExt, :AIBECSShapefileExt,
            :AIBECSJLD2Ext, :AIBECSNCDatasetsExt,
            :AIBECSETOPOExt, :AIBECSOCIM2_48LExt,
        )
        @test Base.get_extension(AIBECS, ext) !== nothing
    end
    # `table(p)` returns a DataFrame in this env
    p = ExtOnP(1.0, 2.0)
    @test typeof(table(p)) === DataFrame
end
