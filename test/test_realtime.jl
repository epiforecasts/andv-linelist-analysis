## Numerical and AD tests for the cluster-completeness integral F_cluster.
##
## The default implementation is a fixed-node Gauss-Hermite quadrature
## (`F_cluster`). `F_cluster_quadrature` is an Integrals.jl/HCubatureJL
## reference solution to ~1e-9. We pin the two against each other and check
## that Enzyme via DifferentiationInterface produces gradients matching
## ForwardDiff to within fixed-node accuracy.

@testset "F_cluster: value agrees with adaptive quadrature reference" begin
    # A spread of (t, μ_inc, σ_inc, μ_δ, σ_δ) covering the regime the model
    # operates in: t between 10 and 80 days, μ_inc ≈ log 20, σ_inc small,
    # μ_δ near zero, σ_δ around 1.
    cases = [
        (15.0, 3.0, 0.4, 0.0,  1.0),
        (30.0, 3.0, 0.5, 0.0,  1.0),
        (50.0, 3.0, 0.5, 1.0,  2.0),
        (80.0, 3.0, 0.6, -2.0, 1.5),
    ]
    for c in cases
        fast = F_cluster(c...)
        ref  = F_cluster_quadrature(c...)
        # Fixed-node Gauss-Hermite at K=40 is accurate to better than 1e-3
        # in this regime; tighter is not needed for likelihood evaluation
        # since the log-likelihood is dominated by other contributions.
        @test isapprox(fast, ref; atol = 1e-3, rtol = 1e-3)
        # Both should be valid probabilities.
        @test 0 <= fast <= 1
    end
end

@testset "F_cluster: monotone in t and limits" begin
    p = (3.0, 0.5, 0.0, 1.0)
    @test F_cluster(0.0,   p...) == 0
    @test F_cluster(-1.0,  p...) == 0
    @test F_cluster(1.0,   p...) < F_cluster(10.0,  p...)
    @test F_cluster(10.0,  p...) < F_cluster(80.0,  p...)
    # Large t -> 1.
    @test F_cluster(300.0, p...) > 0.999
end

@testset "F_cluster: Enzyme gradient via DifferentiationInterface" begin
    θ = [30.0, 3.0, 0.5, 0.0, 1.0]

    # ForwardDiff is the trustworthy reference for a smooth, low-D integrand.
    fwd_val, fwd_grad = value_and_gradient(F_cluster_vec, AutoForwardDiff(), θ)

    # Enzyme reverse mode via DI — this is the deployment target inside the
    # Turing model under AutoEnzyme NUTS.
    enz_backend = AutoEnzyme(; mode = Enzyme.set_runtime_activity(Enzyme.Reverse))
    enz_val, enz_grad = value_and_gradient(F_cluster_vec, enz_backend, θ)

    @test isapprox(fwd_val, enz_val; atol = 1e-12, rtol = 1e-10)
    @test isapprox(fwd_grad, enz_grad; atol = 1e-10, rtol = 1e-8)
    # Sanity: t-derivative must be ≥ 0 (CDF is monotone in its argument).
    @test enz_grad[1] >= 0
end

@testset "F_cluster_quadrature: Enzyme gradient (reference path)" begin
    # The adaptive quadrature path is the one we are NOT putting in the
    # Turing model — but we still want to know if Enzyme can handle it, so
    # this test documents the status. Marked `broken` if Integrals.jl +
    # HCubatureJL is not Enzyme-compatible at this version pin.
    θ = [30.0, 3.0, 0.5, 0.0, 1.0]
    enz_backend = AutoEnzyme(; mode = Enzyme.set_runtime_activity(Enzyme.Reverse))
    try
        _, g = value_and_gradient(F_cluster_quadrature_vec, enz_backend, θ)
        # If it works, check shape and finiteness.
        @test length(g) == 5
        @test all(isfinite, g)
    catch err
        @test_broken false  # record the failure path
        @info "Integrals.jl + Enzyme path not currently working" exception = err
    end
end
