@testset "Particles" begin
    @testset "scalar settling velocity" begin
        T = transportoperator(grd, 100.0)
        @test T isa SparseMatrixCSC
        @test norm(v) / norm(T' * v) > ustrip(upreferred(1u"Myr"))
        T2 = iiptransportoperator(grd; w = p -> (z -> 100.0))
        update_coefficients!(T2, nothing, nothing, nothing)
        @test dropzeros(T2.A) ≈ T
    end



    @testset "settling velocity as a function" begin
        w(z) = upreferred(100u"m/d" + z * u"m" / 100u"m" * 10u"m/d")
        @testset "no sed remin" begin
            T = transportoperator(grd, w)
            @test T isa SparseMatrixCSC
            @test norm(v) / norm(T' * v) > ustrip(upreferred(1u"Myr"))
        end

        @testset "fractional sed remin" begin
            T = transportoperator(grd, w, fsedremin=0.5)
            @test T isa SparseMatrixCSC
            @test norm(v) / norm(T' * v) < ustrip(upreferred(1u"yr"))
        end
    end

    @testset "Inplace transport operators" begin

    end
end
