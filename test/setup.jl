
# Load the circulation and grid
const wet3d, grd, T_Circulation = Circulation.load()

# Define useful constants and arrays
const iwet = indices_of_wet_boxes(wet3d)
const nb = number_of_wet_boxes(wet3d)
const v = vector_of_volumes(wet3d, grd)
const z = vector_of_depths(wet3d, grd)
const ztop = vector_of_top_depths(wet3d, grd)
# And matrices
const DIV = buildDIV(wet3d, iwet, grd)
const Iabove = buildIabove(wet3d, iwet)

@testset "Circulation and grid" begin
    @testset "wet3d" begin
        @test wet3d isa BitArray{3}
    end
    @testset "grd" begin
        @test grd isa Dict
        grd_keys = ["DYT3d", "DZT3d", "dxt"  , "ZT3d" , "DXT3d", "dzt"  , "dyt"  , "ZW3d" , "xt"   , "zt"   , "yt"]
        @testset "has field $k" for k in grd_keys
            @test haskey(grd, k)
        end
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
    # TODO add tests for consistency of wet3d, grid, and T
end

@testset "Constants, vectors, and matrices types" begin
    # Test their types
    @testset "iwet (indices of wet boxes)" begin
        @test iwet isa Vector{<:Int}
        @test length(iwet) == nb
    end
    @testset "nb (number of wet boxes)" begin
        @test nb isa Int
        @test nb == sum(vec(wet3d))
    end
    @testset "v (vector of volumes)" begin
        @test v isa Vector{Float64}
        @test length(v) == nb
        @test all(v .≥ 0)
        @test v == array_of_volumes(grd)[iwet]
    end
    @testset "z (vector of depths of center of boxes)" begin
        @test z isa Vector{Float64}
        @test length(z) == nb
    end
    @test ztop isa Vector{Float64}
    @test DIV isa SparseMatrixCSC
    @test Iabove isa SparseMatrixCSC
end
