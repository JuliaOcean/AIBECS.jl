@info "Running test/AO.jl — AWESOME OCIM download"
@testset begin
    "AWESOME OCIM"
    AO_path = AO.download_and_unpack()
    @test isdir(AO_path)
end
