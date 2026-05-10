@info "Running test/derivatives.jl — ForwardDiff/DualNumbers derivatives" Circulation
# For the conversion to go through ForwardDiff, I need
#Base.convert(FDT::Type{ForwardDiff.Dual{T,V,N}}, p::AIBECS.Parameters) where {T,V,N} = AIBECS.Parameters(convert(Vector{FDT}, vec(p))...)


@testset "Derivatives" begin
    nt = length(T_all)
    n = nt * nb
    @unpack xgeo = p
    x = xgeo * ones(n)
    testp = AIBECS.reconstruct(TestParameters, 2vec(p)) # other value for testing p
    @testset "∇ₓF" begin
        @test ForwardDiff.jacobian(x -> fun(x, p), x) ≈ fun.jac(x, p) rtol = 1.0e-14
        @test ForwardDiff.jacobian(x -> fun(x, p), 2x) ≈ fun.jac(2x, p) rtol = 1.0e-14
        @test ForwardDiff.jacobian(x -> fun(x, p), -x) ≈ fun.jac(-x, p) rtol = 1.0e-14
        @test ForwardDiff.jacobian(x -> fun(x, testp), x) ≈ fun.jac(x, testp) rtol = 1.0e-14
    end
    @testset "∇ₓf" begin
        @test ForwardDiff.jacobian(x -> [f(x, p)], x) ≈ ∇ₓf(x, p) rtol = 1.0e-14
        @test ForwardDiff.jacobian(x -> [f(x, testp)], x) ≈ ∇ₓf(x, testp) rtol = 1.0e-14
        @test ForwardDiff.jacobian(x -> [f(x, p)], 2x) ≈ ∇ₓf(2x, p) rtol = 1.0e-14
    end
    @testset "Analytical ∂Gs" begin
        for xtest in (x, 2x, xgeo * ones(n) .* exp.(cos.(collect(1:n)) / 10))
            J_ad = fun.jac(xtest, p)
            J_an = fun_analytical.jac(xtest, p)
            @test J_an ≈ J_ad rtol = 1.0e-14
            @test nnz(J_an) == nnz(J_ad)
        end
        @test fun_analytical.jac(x, testp) ≈ fun.jac(x, testp) rtol = 1.0e-14
        # Public helper agrees with ForwardDiff at machine precision
        @test check_∂Gs(sms_all, ∂sms_all, x, p, nb).maxabs < 1.0e-12
        # Negative test: a deliberately-wrong ∂Gs (swap two row entries) should disagree
        wrong_∂sms = (
            (∂sms_DIP_∂DOP, ∂sms_DIP_∂DIP, ∂sms_DIP_∂POP),
            ∂sms_all[2],
            ∂sms_all[3],
        )
        @test check_∂Gs(sms_all, wrong_∂sms, x, p, nb).maxabs > 1.0e-10
    end
end
