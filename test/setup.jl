
# Load the circulation and grid
wet3D, grid, T_Circulation_unit = Circulation.load()
T_Circulation = ustrip.(T_Circulation_unit) # strip units for now

# Define useful constants and arrays
iwet = indices_of_wet_boxes(wet3D)
nb = number_of_wet_boxes(wet3D)
v = ustrip.(vector_of_volumes(wet3D, grid)) # strip units for now
z = ustrip.(vector_of_depths(wet3D, grid)) # strip units for now
ztop = ustrip.(vector_of_top_depths(wet3D, grid)) # strip units for now
# And matrices
DIV = buildDIV(wet3D, iwet, grid)
Iabove = buildIabove(wet3D, iwet)

@testset "Circulation and grid" begin
    @testset "wet3D" begin
        @test wet3D isa BitArray{3}
    end
    @testset "grid" begin
        @test grid isa OceanGrid
    end
    @testset "Transport matrix" begin
        T = T_Circulation
        @test T isa SparseMatrixCSC
        @testset "T 1 = 0 (T is divergence free)" begin
            @test norm(ones(nb)) / norm(T * ones(nb)) > ustrip(upreferred(1u"Myr"))
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
        @test nb == sum(vec(wet3D))
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
