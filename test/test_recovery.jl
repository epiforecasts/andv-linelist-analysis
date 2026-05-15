## Simulation-based recovery tests for the submodels of `TransmissionLinelist.jl`.
##
## Each test simulates data from a submodel's own generative assumptions,
## fits the submodel (wrapped in a small Turing model that adds the
## observation likelihood via `@addlogprob!`), and checks that posterior
## medians recover the truth within a justified tolerance.
##
## Sample sizes and chain lengths are intentionally small to keep the
## test suite quick; tolerances are set wide enough to absorb Monte-Carlo
## noise at those settings.

using Random: Random, MersenneTwister
using Distributions: LogNormal, Normal, NegativeBinomial,
                     logpdf, pdf, cdf
using Statistics: median, std, quantile, mean
using Turing: Turing, @model, NUTS, sample, to_submodel
using Integrals: IntegralProblem, QuadGKJL, solve

# Helper: vector of pooled samples for a scalar parameter from a
# Turing/FlexiChains VNChain.
_chain_vec(chn, sym) = vec(collect(chn[sym]))

@testset "Incubation submodel: μ_inc, σ_inc recover from LogNormal draws" begin
    Random.seed!(2025)
    μ_true, σ_true = 3.0, 0.45
    N = 200
    y = rand(LogNormal(μ_true, σ_true), N)

    @model function wrap_inc(y)
        inc ~ to_submodel(TransmissionLinelist.incubation_model(), false)
        Turing.@addlogprob! sum(logpdf.(inc.dist, y))
    end

    chn = sample(wrap_inc(y), NUTS(), 500; progress = false)
    μ_post = _chain_vec(chn, :μ_inc)
    σ_post = _chain_vec(chn, :σ_inc)

    # With N=200 LogNormal draws, the analytic posterior SD on μ is
    # ≈ σ_true / √N ≈ 0.032. A 4-SD posterior band (≈ ±0.13) is a
    # defensible tolerance that survives MCMC noise from 500 NUTS draws.
    @test abs(median(μ_post) - μ_true) < 4 * std(μ_post)
    @test abs(median(σ_post) - σ_true) < 4 * std(σ_post)
    @test abs(median(μ_post) - μ_true) < 0.15
    @test abs(median(σ_post) - σ_true) < 0.15
end

@testset "Transmission-δ submodel: μ_δ, σ_δ recover from Normal draws" begin
    Random.seed!(2026)
    μ_true, σ_true = 1.5, 1.2
    N = 200
    y = rand(Normal(μ_true, σ_true), N)

    @model function wrap_delta(y)
        del ~
        to_submodel(TransmissionLinelist.transmission_delta_model(), false)
        Turing.@addlogprob! sum(logpdf.(del.dist, y))
    end

    chn = sample(wrap_delta(y), NUTS(), 500; progress = false)
    μ_post = _chain_vec(chn, :μ_δ)
    σ_post = _chain_vec(chn, :σ_δ)

    # Posterior SD on μ_δ ≈ σ_true / √N ≈ 0.085. A 4-SD band leaves
    # plenty of margin for MCMC noise at 500 draws.
    @test abs(median(μ_post) - μ_true) < 4 * std(μ_post)
    @test abs(median(σ_post) - σ_true) < 4 * std(σ_post)
    @test abs(median(μ_post) - μ_true) < 0.4
    @test abs(median(σ_post) - σ_true) < 0.3
end

@testset "Random-walk R(t) submodel: posterior bands cover truth" begin
    Random.seed!(2027)
    n_knots = 8
    # Smooth truth: gentle decline in log-R over the knots.
    σ_rw_true = 0.15
    log_R_true = cumsum(vcat(log(1.4),
        σ_rw_true .* randn(n_knots - 1) .- 0.05))
    R_true = exp.(log_R_true)

    # Per-knot NB offspring counts: many cases per bin so the bin-level
    # mean concentrates near the true R, even at moderate dispersion.
    k_nb = 5.0
    cases_per_bin = 60
    Zobs = [rand(NegativeBinomial(k_nb,
                k_nb / (k_nb + R_true[b])))
            for b in 1:n_knots, _ in 1:cases_per_bin]
    counts = vec(sum(Zobs, dims = 2))   # length n_knots
    n_per_bin = fill(cases_per_bin, n_knots)

    @model function wrap_rt(counts, n_per_bin, k_nb)
        rw ~ to_submodel(
            TransmissionLinelist.random_walk_rt_model(length(counts)), false)
        for b in eachindex(counts)
            R = exp(clamp(rw.log_R[b], -50.0, 50.0))
            # Sum of n iid NB(k, p) is NB(n*k, p).
            counts[b] ~ NegativeBinomial(n_per_bin[b] * k_nb,
                k_nb / (k_nb + R))
        end
    end

    chn = sample(wrap_rt(counts, n_per_bin, k_nb), NUTS(), 500;
        progress = false)

    # Reconstruct log_R draws from σ_rw, log_R_init, ε.
    σ_rw = _chain_vec(chn, :σ_rw)
    logR0 = _chain_vec(chn, :log_R_init)
    ε_stack = chn[:ε, stack = true]   # iters × chains × (n_knots - 1)
    n_draws = length(σ_rw)
    log_R_draws = Matrix{Float64}(undef, n_draws, n_knots)
    for d in 1:n_draws
        # Recover (iter, chain) from linear draw index `d`.
        n_iters = size(ε_stack, 1)
        it = ((d - 1) % n_iters) + 1
        ch = ((d - 1) ÷ n_iters) + 1
        ε_d = [ε_stack[it, ch, j] for j in 1:(n_knots - 1)]
        log_R_draws[d, :] = vcat(logR0[d],
            logR0[d] .+ accumulate(+, σ_rw[d] .* ε_d))
    end

    covered = 0
    for b in 1:n_knots
        lo = quantile(log_R_draws[:, b], 0.025)
        hi = quantile(log_R_draws[:, b], 0.975)
        lo <= log_R_true[b] <= hi && (covered += 1)
    end
    # 95% bands on 8 bins from one simulation: cover ≥ 80% (≥ 6/8).
    @test covered >= ceil(Int, 0.8 * n_knots)
end

@testset "ConvolvedDelays cdf: matches QuadGKJL high-precision reference" begin
    inc = LogNormal(3.0, 0.5)
    del = Normal(0.0, 1.0)
    t = 30.0

    # Adaptive QuadGK reference (gold-standard for the value, but not
    # used in production since it doesn't survive Mooncake reverse-mode
    # AD). Same integrand layout as the production GaussLegendre rule
    # but for a scalar `t`.
    f_scalar(δ, p) = (t - δ) > 0 ?
                     cdf(p.inc, t - δ) * pdf(p.del, δ) : 0.0
    prob = IntegralProblem(f_scalar, (-30.0, 30.0),
        (; inc = inc, del = del))
    v_ref = solve(prob, QuadGKJL(); reltol = 1e-12, abstol = 1e-14).u
    v_gl = cdf(ConvolvedDelays(inc, del), [t])[1]
    # Production uses GaussLegendre(n=80), whose absolute error on this
    # integrand is ~1e-6 at t=30. Use a 5e-6 absolute tolerance — tight
    # enough to detect any regression in node count or bounds, loose
    # enough to absorb the fixed-rule's intrinsic truncation error.
    @test isapprox(v_gl, v_ref; atol = 5e-6, rtol = 5e-6)
end
