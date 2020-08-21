

@testset "Aeolian sources" begin
    s_A_2D = AeolianSources.load()
    @testset "$k" for k in keys(s_A_2D)
        v = s_A_2D[k]
        v isa Vector && continue # skip if v is lat/lon
        # Take annual mean
        v_annual = permutedims(dropdims(mean(v, dims=3), dims=3), (2,1))
        # Regrid to OCIM2 grid
        v_regridded = regrid(v_annual, s_A_2D[:lat], s_A_2D[:lon], grd)
        # Paint the top layer
        v_3D = zeros(size(grd)...)
        v_3D[:,:,1] .= ustrip.(upreferred.(v_regridded * u"kg/m^2/s" / grd.δdepth[1]))
        v_vec = v_3D[iwet]
        @test v_vec isa Vector
        @test size(v_vec) == size(iwet)
    end
end

@testset "River sources" begin
    rivers = Rivers.load()
    s = regrid(rivers, grd)
    @test s isa Vector
    @test size(v) == size(iwet)
end

@testset "Groundwater sources" begin
    gws = GroundWaters.load()
    s = regrid(gws, grd)
    @test s isa Vector
    @test size(v) == size(iwet)
end


