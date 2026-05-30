@info "Running test/setup.jl — circulation, grid, transport matrix" Circulation
# Load the circulation and grid. Every shipped `load()` (toy + OCIM* + OCCA)
# returns a plain `SparseMatrixCSC{Float64}` today — the toy circulations
# strip units inside `CG.T_advection`, the JLD2-backed ones strip on load.
grd, T_Circulation = Circulation.load()

# Define useful constants and arrays
iwet = indices_of_wet_boxes(grd)
nb = number_of_wet_boxes(grd)
v = volumevec(grd)
z = depthvec(grd)
ztop = topdepthvec(grd) # strip units for now
# And matrices
DIV = DIVO(grd)
Iabove = buildIabove(grd)

@testset "Circulation and grid" begin
    @testset "wet3D" begin
        wet3D = iswet(grd)
        @test wet3D isa BitArray{3}
        @test wet3D == grd.wet3D
    end
    @testset "grid" begin
        @test grd isa OceanGrid
    end
    @testset "Transport matrix" begin
        T = T_Circulation
        @test T isa SparseMatrixCSC
        println("Divergence timescale: ", norm(ones(nb)) / norm(T * ones(nb)) * u"s" |> u"Myr")
        println("Mass conservation timescale: ", norm(v) / norm(T' * v) * u"s" |> u"Myr")
        @testset "T 1 = 0 (T is divergence free)" begin
            # OCCA's `conservemass = true` correction enforces vᵀT = 0 by
            # subtracting a diagonal that perturbs column sums, so T·1 ≠ 0
            # by construction — see ext/AIBECSJLD2Ext.jl::OCCA.load.
            if nameof(Circulation) === :OCCA
                @test_broken norm(ones(nb)) / norm(T * ones(nb)) > ustrip(upreferred(1u"Myr"))
            else
                @test norm(ones(nb)) / norm(T * ones(nb)) > ustrip(upreferred(1u"Myr"))
            end
        end
        @testset "vᵀ T = 0 (T conserves mass)" begin
            @test norm(v) / norm(T' * v) > ustrip(upreferred(1u"Myr"))
        end
    end
    # TODO add tests for consistency of wet3D, grid, and T
end

@testset "Constants, vectors, and matrices types" begin
    # Test their types
    @testset "iwet (indices of wet boxes)" begin
        @test iwet isa Vector{<:Int}
        @test length(iwet) == nb
    end
    @testset "nb (number of wet boxes)" begin
        @test nb isa Int
        @test nb == sum(vec(grd.wet3D))
    end
    @testset "v (vector of volumes)" begin
        @test v isa Vector{Float64}
        @test length(v) == nb
        @test all(v .≥ 0)
        @test v * u"m^3" == array_of_volumes(grd)[iwet]
    end
    @testset "z (vector of depths of center of boxes)" begin
        @test z isa Vector{<:Real}
        @test length(z) == nb
    end
    @test ztop isa Vector{<:Real}
    @test DIV isa SparseMatrixCSC
    @test Iabove isa SparseMatrixCSC
end
