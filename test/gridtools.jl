@testset "GridTools" begin
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n)
    @testset "convert to 3D array" for xᵢ in state_to_tracers(x, nb, nt)
        xᵢ3D = rearrange_into_3Darray(xᵢ, grd)
        @test size(xᵢ3D) == size(grd)
    end
    @testset "convert to 3D array" for xᵢ in unpack_tracers(x, grd)
        xᵢ3D = rearrange_into_3Darray(xᵢ, grd)
        @test size(xᵢ3D) == size(grd)
    end
end
