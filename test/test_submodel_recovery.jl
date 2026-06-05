## Additional simulation-based recovery tests using Turing's
## `fix` / `condition` machinery and `rand(model)` for data simulation.
##
## Existing recovery tests in `test_recovery.jl` cover the marginal
## delay submodels (`incubation_model`, `transmission_delta_model`) and
## the random walk (`random_walk_rt_model`). Here we add:
##
## - `nb_dispersion_model`: recover dispersion `k` from simulated NB
##   offspring counts.
## - `case_model`: recover `k`, `log_R_init`, `σ_rw` jointly from a
##   simulated cohort of cases with known onset times.
## - `random_walk_rt_model` (deepened): explicit CrI coverage on
##   `σ_rw` and `log_R_init`, separate from the existing per-knot
##   coverage check.
##
## Submodels skipped (with reason):
##
## - `truncation_model`: no parameters of its own; contributions are
##   `@addlogprob!` terms only. The parameters that would be recovered
##   are the parent's delay distributions, already tested elsewhere.
## - `latent_times_model`: no population parameters of its own; it
##   draws per-case latents from Uniform priors and adds marginal
##   incubation/transmission log-densities given fixed parent
##   `inc_dist`/`δ_dist`. Sim-then-recover would have to happen at
##   the joint level.
## - `delays_only_model`, `joint_model`: composite. `test_recovery.jl`
##   exercises these via convergence and predictor smoke tests; full
##   joint sim-and-recover would balloon test runtime.

using Random: Random
using Distributions: LogNormal, Normal, NegativeBinomial, Uniform,
                     logpdf, truncated
using Statistics: median, std, quantile, mean
using Turing: Turing, @model, NUTS, sample, to_submodel, DynamicPPL

# `_chain_vec` is defined in `test_recovery.jl`, included earlier in
# `runtests.jl`, so it is in scope here.

@testset "nb_dispersion_model: k recovers from simulated NB counts" begin
    Random.seed!(2030)
    k_true = 2.5
    phi_inv_sqrt_true = 1.0 / sqrt(k_true)
    R_true = 1.2
    N = 400

    # Step 1: wrap `nb_dispersion_model` in a small Turing model that
    # draws NB counts with known `R`. Fix `phi_inv_sqrt` to truth so
    # `rand(model)` simulates `Z` from the fixed-dispersion likelihood.
    @model function sim_nb(N, R)
        disp ~ to_submodel(TransmissionLinelist.nb_dispersion_model(), false)
        Z = Vector{Int}(undef, N)
        for i in 1:N
            Z[i] ~ NegativeBinomial(disp.k, disp.k / (disp.k + R))
        end
        return Z
    end
    sim_model = DynamicPPL.fix(sim_nb(N, R_true),
        (; phi_inv_sqrt = phi_inv_sqrt_true))
    sim = rand(sim_model)
    Z_sim = [sim[k] for k in keys(sim)]

    # Step 2: build the same model with concrete observed `Z` (which
    # `@model` treats as data) and sample. Equivalent to using
    # `condition` to pin `Z` to the simulated draws.
    @model function fit_nb(Z, R)
        disp ~ to_submodel(TransmissionLinelist.nb_dispersion_model(), false)
        for i in eachindex(Z)
            Z[i] ~ NegativeBinomial(disp.k, disp.k / (disp.k + R))
        end
    end
    chn = sample(fit_nb(Z_sim, R_true), NUTS(), 500; progress = false)
    k_post = _chain_vec(chn, :k)

    lo, hi = quantile(k_post, [0.025, 0.975])
    @test lo <= k_true <= hi
    @test abs(median(k_post) - k_true) < 4 * std(k_post)
end

@testset "case_model: log_R_init, σ_rw, k recover from simulated cases" begin
    Random.seed!(2031)
    n_knots = 6
    # Knot positions, 7 days apart.
    edges = collect(0.0:7.0:(7.0 * (n_knots - 1)))
    # Cases spread across the window.
    N = 240
    T_onset = collect(range(0.5, edges[end] - 0.5; length = N))
    # Retrospective mode: empty `p`.
    p = Float64[]

    # Truth for the rt + dispersion submodels.
    σ_rw_true = 0.12
    log_R_init_true = log(1.3)
    ε_true = randn(n_knots - 1)
    phi_inv_sqrt_true = 1.0 / sqrt(4.0)   # k_true = 4

    truth = (; σ_rw = σ_rw_true,
        log_R_init = log_R_init_true,
        ε = ε_true,
        phi_inv_sqrt = phi_inv_sqrt_true)

    # Step 1: build `case_model` with `Z = Vector{Missing}` so `Z` is
    # treated as latent. Fix the priors to truth, `rand` to simulate.
    sim_case = TransmissionLinelist.case_model(
        Vector{Missing}(missing, N), edges, T_onset, p)
    sim_fixed = DynamicPPL.fix(sim_case, truth)
    sim = rand(sim_fixed)
    Z_sim = [Int(sim[k]) for k in keys(sim)]

    # Step 2: build a fresh `case_model` with the simulated counts as
    # concrete observed data and sample. `case_model` treats `Z` as
    # observed when given a concrete integer vector.
    fit_case = TransmissionLinelist.case_model(Z_sim, edges, T_onset, p)
    chn = TransmissionLinelist.sample_fit(fit_case;
        samples = 500, chains = 2, seed = 20260518, progress = false)

    σ_rw_post = _chain_vec(chn, :σ_rw)
    logR0_post = _chain_vec(chn, :log_R_init)
    k_post = _chain_vec(chn, :k)

    for (post,
        truth_val) in [(σ_rw_post, σ_rw_true),
        (logR0_post, log_R_init_true),
        (k_post, 4.0)]
        lo, hi = quantile(post, [0.025, 0.975])
        @test lo <= truth_val <= hi
    end
end

@testset "random_walk_rt_model: σ_rw and log_R_init CrI cover truth" begin
    # Deepens the existing per-knot coverage test by checking that the
    # population-level parameters (`σ_rw`, `log_R_init`) themselves
    # sit inside their 95% posterior CrIs.
    Random.seed!(2032)
    n_knots = 10
    σ_rw_true = 0.15
    log_R_init_true = log(1.4)
    ε_true = randn(n_knots - 1)
    log_R_true = vcat(log_R_init_true,
        log_R_init_true .+ accumulate(+, σ_rw_true .* ε_true))

    # Step 1: simulate per-knot NB counts using the truth.
    k_nb = 5.0
    cases_per_bin = 80
    counts = [sum(rand(
                  NegativeBinomial(k_nb,
                      k_nb / (k_nb + exp(log_R_true[b]))),
                  cases_per_bin)) for b in 1:n_knots]
    n_per_bin = fill(cases_per_bin, n_knots)

    # Step 2: fit the rt submodel to the simulated counts. Same wrap
    # as the existing per-knot test, but we only inspect scalar pop
    # params here.
    @model function wrap_rt(counts, n_per_bin, k_nb)
        rw ~ to_submodel(
            TransmissionLinelist.random_walk_rt_model(length(counts)), false)
        for b in eachindex(counts)
            R = exp(clamp(rw.log_R[b], -50.0, 50.0))
            counts[b] ~ NegativeBinomial(n_per_bin[b] * k_nb,
                k_nb / (k_nb + R))
        end
    end
    chn = sample(wrap_rt(counts, n_per_bin, k_nb), NUTS(), 500;
        progress = false)
    σ_rw_post = _chain_vec(chn, :σ_rw)
    logR0_post = _chain_vec(chn, :log_R_init)

    lo, hi = quantile(σ_rw_post, [0.025, 0.975])
    @test lo <= σ_rw_true <= hi
    lo, hi = quantile(logR0_post, [0.025, 0.975])
    @test lo <= log_R_init_true <= hi
end
