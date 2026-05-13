## Tests for F_offspring: numerical value tracks a Monte Carlo
## reference, and the function is differentiable under Mooncake reverse
## mode via DifferentiationInterface.jl (the AD backend joint_model uses).
##
## The Mooncake test exercises the multi-element-vector call shape that
## joint_model itself uses, so coverage matches production.

using Random: MersenneTwister
using Distributions: LogNormal, Normal

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
        v_qr = F_offspring([t], LogNormal(μ_inc, σ_inc), Normal(μ_δ, σ_δ))[1]
        v_mc = _f_offspring_mc(t, μ_inc, σ_inc, μ_δ, σ_δ)
        @test isapprox(v_qr, v_mc; atol = 5e-3)
        @test 0 <= v_qr <= 1
    end
end

@testset "F_offspring: monotone in t and limits" begin
    inc = LogNormal(3.0, 0.5)
    del = Normal(0.0, 1.0)
    F(t) = F_offspring([t], inc, del)[1]
    @test F(0.0)   <= 1e-6
    @test F(-50.0) <= 1e-10
    @test F(1.0)   < F(10.0)
    @test F(10.0)  < F(80.0)
    @test F(300.0) > 0.999
end

@testset "F_offspring: Mooncake reverse via DifferentiationInterface" begin
    ts = collect(10.0:5.0:55.0)   # 10 points — multi-element, matches joint_model use
    loss(θ) = sum(F_offspring(ts, LogNormal(θ[1], θ[2]), Normal(θ[3], θ[4])))
    θ = [3.0, 0.5, 0.0, 1.0]

    fwd_val, fwd_grad = value_and_gradient(loss, AutoForwardDiff(), θ)

    mc_backend = AutoMooncake(; config = Mooncake.Config())
    mc_val, mc_grad = value_and_gradient(loss, mc_backend, θ)

    @test isapprox(fwd_val, mc_val; atol = 1e-10, rtol = 1e-10)
    @test isapprox(fwd_grad, mc_grad; atol = 1e-8, rtol = 1e-6)
end
