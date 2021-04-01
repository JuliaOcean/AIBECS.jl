
# Load the circulation and grid
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
    println("    Circulation and grid")
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
        @test T isa DiffEqArrayOperator
        @test T.A isa SparseMatrixCSC
        println("      Divergence timescale: ", norm(ones(nb)) / norm(T * ones(nb)) * u"s" |> u"Myr")
        println("      Mass conservation timescale: ", norm(v) / norm(T.A' * v) * u"s" |> u"Myr") # TODO adjoint(Op)
    end
    # TODO add tests for consistency of wet3D, grid, and T
end

@testset "Constants, vectors, and matrices types" begin
    println("    Constants, vectors, and matrices types")
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

