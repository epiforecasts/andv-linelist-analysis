## Tests for F_cluster: it is differentiable under Mooncake reverse
## mode via DifferentiationInterface.jl (the AD backend joint_model
## uses), and its numerical value tracks a Monte Carlo reference for
## the underlying probability.

using Random: MersenneTwister
using Distributions: LogNormal, Normal

# Brute-force MC reference: draws (Inc(src), δ, Inc(sec)) from the
# joint and counts the fraction with sum ≤ t. The QR default should
# track this within MC noise (1e4 draws gives ~1e-2 SE; we test 1e-2
# tolerance).
function _f_cluster_mc(t, μ_inc, σ_inc, μ_δ, σ_δ; n = 100_000, seed = 42)
    rng = MersenneTwister(seed)
    inc = LogNormal(μ_inc, σ_inc)
    del = Normal(μ_δ, σ_δ)
    hits = 0
    for _ in 1:n
        x1 = rand(rng, inc); y = rand(rng, del); x2 = rand(rng, inc)
        x1 + y + x2 <= t && (hits += 1)
    end
    return hits / n
end

@testset "F_cluster: value tracks a Monte Carlo reference" begin
    cases = [
        (30.0, 3.0, 0.5, 0.0,  1.0),
        (50.0, 3.0, 0.5, 1.0,  2.0),
        (80.0, 3.0, 0.6, -2.0, 1.5),
    ]
    for c in cases
        v_qr = F_cluster(c...)
        v_mc = _f_cluster_mc(c...)
        @test isapprox(v_qr, v_mc; atol = 5e-3)
        @test 0 <= v_qr <= 1
    end
end

@testset "F_cluster: monotone in t and limits" begin
    p = (3.0, 0.5, 0.0, 1.0)
    @test F_cluster(0.0,   p...) <= 1e-10
    @test F_cluster(-1.0,  p...) <= 1e-10
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
    @test mc_grad[1] >= 0   # CDF is monotone in t
end
