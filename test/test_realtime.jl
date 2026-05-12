## Tests for F_cluster: it is differentiable under Mooncake reverse mode
## via DifferentiationInterface.jl (the AD backend joint_model uses) and
## its value is invariant to the choice of Integrals.jl algorithm.

using Integrals: HCubatureJL, QuadratureRule

@testset "F_cluster: value invariant under Integrals.jl algorithm" begin
    cases = [
        (15.0, 3.0, 0.4, 0.0,  1.0),
        (30.0, 3.0, 0.5, 0.0,  1.0),
        (50.0, 3.0, 0.5, 1.0,  2.0),
        (80.0, 3.0, 0.6, -2.0, 1.5),
    ]
    for c in cases
        v_default = F_cluster(c...)
        v_loose   = F_cluster(c...; reltol = 1e-4, abstol = 1e-6)
        @test isapprox(v_default, v_loose; atol = 1e-4, rtol = 1e-4)
        @test 0 <= v_default <= 1
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

@testset "F_cluster: Mooncake reverse via DifferentiationInterface" begin
    θ = [30.0, 3.0, 0.5, 0.0, 1.0]

    fwd_val, fwd_grad = value_and_gradient(F_cluster_vec, AutoForwardDiff(), θ)

    mc_backend = AutoMooncake(; config = Mooncake.Config())
    mc_val, mc_grad = value_and_gradient(F_cluster_vec, mc_backend, θ)

    @test isapprox(fwd_val, mc_val; atol = 1e-10, rtol = 1e-10)
    @test isapprox(fwd_grad, mc_grad; atol = 1e-8, rtol = 1e-6)
    @test mc_grad[1] >= 0  # CDF is monotone in t
end
