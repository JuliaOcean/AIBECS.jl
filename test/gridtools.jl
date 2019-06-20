@testset "GridTools" begin
    nt = length(T_all)
    n = nt * nb
    x₀ = p₀.xgeo * ones(n)
    @testset "convert to 3D array" for xᵢ in state_to_tracers(x₀, nb, nt)
        xᵢ3D = rearrange_into_3Darray(xᵢ, wet3D)
        @test size(xᵢ3D) == size(wet3D)
        @test size(rearrange_into_1Dvector(xᵢ3D, wet3D)) == size(xᵢ)
    end
end
