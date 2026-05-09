@info "Running test/sparse_jacobian.jl — bespoke vs DI+SCT Jacobian" Circulation
using DifferentiationInterface
using SparseConnectivityTracer
using SparseMatrixColorings
using ADTypes
using SparseArrays

# Validates AIBECS's bespoke sparse Jacobian (`prob.f.jac`, assembled in
# `src/multiTracer.jl:local_jacobian`) against a reference Jacobian computed via
# DifferentiationInterface + SparseConnectivityTracer + ForwardDiff.

@testset "Sparse Jacobian: AIBECS bespoke vs DI+SCT" begin
    backend = AutoSparse(
        AutoForwardDiff();
        sparsity_detector = TracerSparsityDetector(),
        coloring_algorithm = GreedyColoringAlgorithm(),
    )

    builders = (
        idealage    = build_idealage_problem,
        radiocarbon = build_radiocarbon_problem,
        po4pop      = build_pmodel_problem,
    )

    for (tname, build) in pairs(builders)
        @testset "$tname" begin
            prob = build(grd, T_Circulation)
            p    = prob.p
            u0   = prob.u0
            F    = x -> prob.f(x, p, 0.0)

            for (xlabel, x) in (("u0", u0), ("perturbed", 1.3 .* u0 .+ 1.0e-3))
                @testset "x = $xlabel" begin
                    prep     = prepare_jacobian(F, backend, x)
                    J_aibecs = dropzeros(prob.f.jac(x, p))
                    J_di     = dropzeros(jacobian(F, prep, backend, x))

                    @test J_aibecs ≈ J_di rtol = 1.0e-12
                    @test nnz(J_aibecs) == nnz(J_di)
                end
            end
        end
    end
end
