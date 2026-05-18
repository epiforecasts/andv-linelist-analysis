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
                     logpdf, pdf, cdf, ccdf
using Dates: Dates, Date, Day
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

@testset "delays_only_model: real-time mode converges at both cut-offs" begin
    # Regression guard for the σ_δ collapse and Mooncake bitcast classes
    # of bug: if cdf(ConvolvedDelays, ·) ever stops being AD-stable
    # under parametric (μ, σ) the real-time fit blows up here.
    ll = TransmissionLinelist.load_linelist()
    t0 = minimum(ll.onset_date) - Day(60)
    for obs_date in [Date("2018-12-31"), Date("2019-01-07")]
        ll_rt = TransmissionLinelist.filter_realtime(ll, obs_date)
        d_rt = TransmissionLinelist.build_data(ll_rt;
            obs_time = obs_date, t0 = t0)
        chn = TransmissionLinelist.sample_fit(
            TransmissionLinelist.delays_only_model(d_rt);
            samples = 200, chains = 2, seed = 20260512, progress = false)
        diag = TransmissionLinelist.diagnostics_table(chn)
        @test diag.rhat_max[1] < 1.1
        @test diag.divergences[1] / (2 * 200) < 0.05
    end
end

@testset "predict_*_outbreak: shapes + strict ≤ natural" begin
    # Smoke test for both real-time predictors. Strict counterfactual
    # (no further transmission) should sit below the natural-chain one
    # in expectation: the per-source probability `q_i − p_i` ≤ `1 − p_i`.
    ll = TransmissionLinelist.load_linelist()
    t0 = minimum(ll.onset_date) - Day(60)
    obs_date = Date("2018-12-31")
    ll_rt = TransmissionLinelist.filter_realtime(ll, obs_date)
    d_rt = TransmissionLinelist.build_data(ll_rt;
        obs_time = obs_date, t0 = t0)
    edges = TransmissionLinelist.prepare_rt_edges(t0; obs_time = obs_date)
    m = TransmissionLinelist.joint_model(d_rt, edges)
    chn = TransmissionLinelist.sample_fit(m;
        samples = 200, chains = 2, seed = 20260512, progress = false)
    post = TransmissionLinelist.summarise(chn)
    strict = TransmissionLinelist.predict_controlled_outbreak(
        m, chn, post, d_rt; obs_time = obs_date, t0 = t0)
    natural = TransmissionLinelist.predict_natural_chain_outbreak(
        m, chn, post, d_rt; obs_time = obs_date, t0 = t0)
    @test length(strict.future_samples) == length(natural.future_samples)
    @test all(>=(0), strict.future_samples)
    @test all(>=(0), natural.future_samples)
    @test TransmissionLinelist.realised_future_count(ll, obs_date) >= 0
    @test mean(strict.future_samples) <= mean(natural.future_samples) + 1
end

@testset "predict_controlled_outbreak: intervention_time variants" begin
    # Scalar Date earlier than the observation cut-off reduces the
    # mean future count vs. the default (`intervention_time = nothing`,
    # equivalent to `obs_date`).
    ll = TransmissionLinelist.load_linelist()
    t0 = minimum(ll.onset_date) - Day(60)
    obs_date = Date("2018-12-31")
    ll_rt = TransmissionLinelist.filter_realtime(ll, obs_date)
    d_rt = TransmissionLinelist.build_data(ll_rt;
        obs_time = obs_date, t0 = t0)
    edges = TransmissionLinelist.prepare_rt_edges(t0; obs_time = obs_date)
    m = TransmissionLinelist.joint_model(d_rt, edges)
    chn = TransmissionLinelist.sample_fit(m;
        samples = 100, chains = 2, seed = 20260512, progress = false)
    post = TransmissionLinelist.summarise(chn)

    earlier = obs_date - Day(7)
    s_scalar = TransmissionLinelist.predict_controlled_outbreak(
        m, chn, post, d_rt;
        obs_time = obs_date, t0 = t0,
        intervention_time = earlier,
        rng = Random.MersenneTwister(42))
    s_obs = TransmissionLinelist.predict_controlled_outbreak(
        m, chn, post, d_rt;
        obs_time = obs_date, t0 = t0,
        rng = Random.MersenneTwister(42))
    @test mean(s_scalar.future_samples) <= mean(s_obs.future_samples)
end

@testset "MC validation: _pipeline_probability matches Monte Carlo" begin
    # Validate `_pipeline_probability` against direct Monte Carlo for
    # the joint event `{δ ≤ Δ_q ∧ δ + Inc > Δ_p}`. Covers the case
    # that regressed under the old `q − p` formula (intervention
    # strictly before the observation horizon, `Δ_q < Δ_p`), the
    # natural-chain limit (`Δ_q = +Inf`, must equal `1 - cdf(
    # ConvolvedDelays(inc, δd), Δ_p)`), and the exclusion sentinel
    # (`Δ_q = -Inf` must return exactly zero).
    rng = Random.MersenneTwister(20260518)
    N_mc = 10_000
    δ_dist = Normal(8.0, 4.0)
    inc_dist = LogNormal(2.0, 0.5)
    δ_draws = rand(rng, δ_dist, N_mc)
    inc_draws = rand(rng, inc_dist, N_mc)
    sum_draws = δ_draws .+ inc_draws

    cases = [
        (Δ_q = 10.0, Δ_p = 10.0, label = "Δq=Δp (regression check)"),
        (Δ_q = 5.0, Δ_p = 10.0, label = "Δq<Δp by 5d (old formula off)"),
        (Δ_q = 0.0, Δ_p = 10.0, label = "Δq=0"),
        (Δ_q = -2.0, Δ_p = 10.0, label = "Δq<0")
    ]
    for c in cases
        joint_emp = mean((δ_draws .<= c.Δ_q) .& (sum_draws .> c.Δ_p))
        mcse = sqrt(joint_emp * (1 - joint_emp) / N_mc)
        analytic = TransmissionLinelist._pipeline_probability(
            inc_dist, δ_dist, c.Δ_q, c.Δ_p)
        @test isapprox(analytic, joint_emp; atol = 5 * mcse + 1e-3)
    end

    # Δ_q = +Inf must reduce to 1 − cdf(ConvolvedDelays, Δ_p).
    for Δ_p in (10.0, 30.0)
        conv = TransmissionLinelist.ConvolvedDelays(inc_dist, δ_dist)
        natural = TransmissionLinelist._pipeline_probability(
            inc_dist, δ_dist, Inf, Δ_p)
        @test isapprox(natural, 1 - cdf(conv, Δ_p); atol = 1e-6)
    end

    # Δ_q = -Inf must return exactly zero.
    @test TransmissionLinelist._pipeline_probability(
        inc_dist, δ_dist, -Inf, 10.0) == 0.0
end
