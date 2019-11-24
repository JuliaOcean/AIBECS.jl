@testset "GridTools" begin
    nt = length(T_all)
    n = nt * nb
    x₀ = p₀.xgeo * ones(n)
    @testset "convert to 3D array" for xᵢ in state_to_tracers(x₀, nb, nt)
        xᵢ3D = rearrange_into_3Darray(xᵢ, grid)
        @test size(xᵢ3D) == size(grid)
    end
    @testset "convert to 3D array" for xᵢ in state_to_tracers(x₀, grid)
        xᵢ3D = rearrange_into_3Darray(xᵢ, grid)
        @test size(xᵢ3D) == size(grid)
    end
end
