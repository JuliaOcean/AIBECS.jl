
# Load the circulation and grid
grid, T_Circulation_unit = Circulation.load()
T_Circulation = ustrip.(T_Circulation_unit) # strip units for now

# Define useful constants and arrays
iwet = indices_of_wet_boxes(grid)
nb = number_of_wet_boxes(grid)
v = ustrip.(vector_of_volumes(grid)) # strip units for now
z = ustrip.(vector_of_depths(grid)) # strip units for now
ztop = ustrip.(vector_of_top_depths(grid)) # strip units for now
# And matrices
DIV = buildDIV(grid)
Iabove = buildIabove(grid)

@testset "Circulation and grid" begin
    @testset "wet3D" begin
        wet3D = iswet(grid)
        @test wet3D isa BitArray{3}
        @test wet3D == grid.wet3D
    end
    @testset "grid" begin
        @test grid isa OceanGrid
    end
    @testset "Transport matrix" begin
        T = T_Circulation
        @test T isa SparseMatrixCSC
        println("Divergence timescale: ", norm(ones(nb)) / norm(T * ones(nb)) * u"s" |> u"Myr")
        println("Mass conservation timescale: ", norm(v) / norm(T' * v) * u"s" |> u"Myr")
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
        @test nb == sum(vec(grid.wet3D))
    end
    @testset "v (vector of volumes)" begin
        @test v isa Vector{Float64}
        @test length(v) == nb
        @test all(v .≥ 0)
        @test v * u"m^3" == array_of_volumes(grid)[iwet]
    end
    @testset "z (vector of depths of center of boxes)" begin
        @test z isa Vector{<:Real}
        @test length(z) == nb
    end
    @test ztop isa Vector{<:Real}
    @test DIV isa SparseMatrixCSC
    @test Iabove isa SparseMatrixCSC
end

