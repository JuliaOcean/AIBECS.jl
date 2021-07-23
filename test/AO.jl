@testset begin "AWESOME OCIM"
    AO.download_and_unpack()
    @test isdir(datadep"AWESOME-OCIM")
end
