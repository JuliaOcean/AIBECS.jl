@testset begin "AWESOME OCIM"
    AO.download_and_unpack()
    @test isdir(joinpath(homedir(), ".julia", "datadeps", "AWESOME-OCIM"))
end
