## Tests for F_offspring: it is differentiable under Mooncake reverse
## mode via DifferentiationInterface.jl (the AD backend joint_model
## uses), and its numerical value tracks a Monte Carlo reference for
## the underlying probability.

using Random: MersenneTwister
using Distributions: LogNormal, Normal

# Brute-force MC reference: draws (δ, Inc(sec)) from the joint and
# counts the fraction with sum ≤ t. The source's own incubation is not
# part of the offspring-completeness probability — it is conditioned on
# inside `joint_model` via the sampled `T_onset[src]` latent. The QR
# default should track this within MC noise (1e5 draws gives ~1e-3 SE;
# we test 5e-3 tolerance).
function _f_offspring_mc(t, μ_inc, σ_inc, μ_δ, σ_δ; n = 100_000, seed = 42)
    rng = MersenneTwister(seed)
    inc = LogNormal(μ_inc, σ_inc)
    del = Normal(μ_δ, σ_δ)
    hits = 0
    for _ in 1:n
        y = rand(rng, del); x2 = rand(rng, inc)
        y + x2 <= t && (hits += 1)
    end
    return hits / n
end

@testset "F_offspring: value tracks a Monte Carlo reference" begin
    cases = [
        (30.0, 3.0, 0.5, 0.0,  1.0),
        (50.0, 3.0, 0.5, 1.0,  2.0),
        (80.0, 3.0, 0.6, -2.0, 1.5),
    ]
    for (t, μ_inc, σ_inc, μ_δ, σ_δ) in cases
        v_qr = F_offspring(t, LogNormal(μ_inc, σ_inc), Normal(μ_δ, σ_δ))
        v_mc = _f_offspring_mc(t, μ_inc, σ_inc, μ_δ, σ_δ)
        @test isapprox(v_qr, v_mc; atol = 5e-3)
        @test 0 <= v_qr <= 1
    end
end

@testset "F_offspring: monotone in t and limits" begin
    inc = LogNormal(3.0, 0.5)
    del = Normal(0.0, 1.0)
    @test F_offspring(0.0,   inc, del) <= 1e-6
    @test F_offspring(-50.0, inc, del) <= 1e-10
    @test F_offspring(1.0,   inc, del) < F_offspring(10.0, inc, del)
    @test F_offspring(10.0,  inc, del) < F_offspring(80.0, inc, del)
    @test F_offspring(300.0, inc, del) > 0.999
end

@testset "F_offspring: Mooncake reverse via DifferentiationInterface" begin
    θ = [30.0, 3.0, 0.5, 0.0, 1.0]

    fwd_val, fwd_grad = value_and_gradient(F_offspring_vec, AutoForwardDiff(), θ)

    mc_backend = AutoMooncake(; config = Mooncake.Config())
    mc_val, mc_grad = value_and_gradient(F_offspring_vec, mc_backend, θ)

    @test isapprox(fwd_val, mc_val; atol = 1e-10, rtol = 1e-10)
    @test isapprox(fwd_grad, mc_grad; atol = 1e-8, rtol = 1e-6)
    @test mc_grad[1] >= 0   # CDF is monotone in t
end
