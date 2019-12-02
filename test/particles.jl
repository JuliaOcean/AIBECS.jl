@testset "Particles" begin
    T = transportoperator(grid, w=100.0)
    @test T isa SparseMatrixCSC
    @testset "vᵀ T = 0 (T conserves mass)" begin
        @test norm(v) / norm(T' * v) > ustrip(upreferred(1u"Myr"))
    end
end
