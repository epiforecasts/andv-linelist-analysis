## Numerical and AD tests for the cluster-completeness integral F_cluster.
##
## Value tests compare the default Gauss-Hermite path to an Integrals.jl
## HCubatureJL reference (without AD). The AD test exercises only the
## production path (Gauss-Hermite under Enzyme reverse mode via
## DifferentiationInterface), since reverse-mode Enzyme does not work
## through `Integrals.solve` at the current SciMLBase/Enzyme pin
## regardless of the underlying algorithm.

@testset "F_cluster: value agrees with Integrals.jl HCubature reference" begin
    cases = [
        (15.0, 3.0, 0.4, 0.0,  1.0),
        (30.0, 3.0, 0.5, 0.0,  1.0),
        (50.0, 3.0, 0.5, 1.0,  2.0),
        (80.0, 3.0, 0.6, -2.0, 1.5),
    ]
    for c in cases
        fast = F_cluster(c...)
        ref  = F_cluster(c...; alg = HCubature())
        # Default Gauss-Hermite is accurate to better than 1e-3 against
        # the adaptive HCubature reference across the model regime.
        @test isapprox(fast, ref; atol = 1e-3, rtol = 1e-3)
        @test 0 <= fast <= 1
    end
end

@testset "F_cluster: monotone in t and limits" begin
    p = (3.0, 0.5, 0.0, 1.0)
    @test F_cluster(0.0,   p...) == 0
    @test F_cluster(-1.0,  p...) == 0
    @test F_cluster(1.0,   p...) < F_cluster(10.0,  p...)
    @test F_cluster(10.0,  p...) < F_cluster(80.0,  p...)
    @test F_cluster(300.0, p...) > 0.999
end

@testset "F_cluster: GaussHermite converges to HCubature as n grows" begin
    c = (50.0, 3.0, 0.5, 1.0, 2.0)
    ref = F_cluster(c...; alg = HCubature())
    err20 = abs(F_cluster(c...; alg = GaussHermite(20)) - ref)
    err40 = abs(F_cluster(c...; alg = GaussHermite(40)) - ref)
    err80 = abs(F_cluster(c...; alg = GaussHermite(80)) - ref)
    @test err40 < err20
    @test err80 < 1e-3
end

@testset "F_cluster (GaussHermite): Enzyme reverse via DifferentiationInterface" begin
    θ = [30.0, 3.0, 0.5, 0.0, 1.0]

    fwd_val, fwd_grad = value_and_gradient(F_cluster_vec, AutoForwardDiff(), θ)

    enz_backend = AutoEnzyme(; mode = Enzyme.set_runtime_activity(Enzyme.Reverse))
    enz_val, enz_grad = value_and_gradient(F_cluster_vec, enz_backend, θ)

    @test isapprox(fwd_val, enz_val; atol = 1e-12, rtol = 1e-10)
    @test isapprox(fwd_grad, enz_grad; atol = 1e-10, rtol = 1e-8)
    @test enz_grad[1] >= 0   # CDF is monotone in t
end
