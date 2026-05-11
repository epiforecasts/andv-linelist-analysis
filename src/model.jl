## Joint model for the Epuyén ANDV outbreak: incubation period, per-pair
## transmission timing relative to source onset, and time-varying reproduction
## number.
##
## Three quantities estimated together from the line list:
##
##   1. Incubation period           — LogNormal(μ_inc, σ_inc), from each case's
##                                    exposure-window-to-onset gap.
##   2. Transmission timing δ       — Normal(μ_δ, σ_δ), the gap between a
##                                    secondary's infection time and its
##                                    source's symptom onset. Identified
##                                    per-pair from the line list.
##   3. Time-varying reproduction   — log R(t) on a weekly random walk, with
##      number                        Negative-Binomial offspring (dispersion k).
##
## Each case has continuous latents: an infection time T_inf and an onset time
## T_onset. Interval-censored onsets and exposure dates are handled by
## Bayesian data augmentation over these latents.
##
## The generation interval and serial interval are derived in post-processing
## from δ and Inc. The per-pair constraint T_inf[secondary] > T_inf[source]
## is enforced via a -Inf reject in the likelihood to ensure GI > 0.

@model function joint_model(d, edges)
    # Population-level parameters
    μ_inc ~ Normal(3.0, 0.5)                       # log-mean Inc (≈ log 20 d)
    σ_inc ~ truncated(Normal(0.0, 0.5); lower = 0) # log-SD Inc
    μ_δ   ~ Normal(0.0, 5.0)                       # population mean transmission timing (d from source onset)
    σ_δ   ~ truncated(Normal(0.0, 1.0); lower = 0) # population SD of transmission timing (d)
    k     ~ truncated(Normal(0.3, 0.5); lower = 0) # NB offspring dispersion (centred low — known super-spreader pathogen)
    σ_rw  ~ truncated(Normal(0.0, 0.5); lower = 0) # log-R RW innovation SD

    # Concrete element type derived from a sampled scalar: stable Float64 on
    # the forward pass and a stable Dual / tracked type under AD. Helps every
    # backend (ForwardDiff, ReverseDiff, Mooncake, Enzyme) — replacing
    # `Vector{Real}` removes a dynamic-dispatch tax inside the inner loop.
    T = typeof(μ_inc)

    # Random walk on log R(t) across time bins. Non-centred parameterisation:
    # sample i.i.d. standard-normal innovations ε, then reconstruct log_R via
    # a cumulative sum scaled by σ_rw. The mapping has unit Jacobian (linear
    # in ε with σ_rw fixed in a step), so the implied prior on the centred
    # log_R vector is unchanged. Removes the σ_rw–log_R funnel that drove
    # ~0.6% of NUTS samples to diverge under the centred form.
    n_bins = length(edges) + 1
    log_R_init ~ Normal(log(1.5), 1.0)
    ε ~ Turing.filldist(Normal(zero(T), one(T)), n_bins - 1)
    log_R := vcat(log_R_init, log_R_init .+ accumulate(+, σ_rw .* ε))

    inc_dist = LogNormal(μ_inc, σ_inc)

    # T_onset is a latent over the recorded onset window (defaults to a
    # one-day window when only a single onset date was recorded).
    T_onset = Vector{T}(undef, d.N)
    for i in 1:d.N
        T_onset[i] ~ Uniform(d.onset_lo_day[i], d.onset_hi_day[i])
    end

    T_inf = Vector{T}(undef, d.N)
    for i in 1:d.N
        if d.source_idx[i] == 0
            # Zoonotic index: free latent T_inf pre-onset.
            T_inf[i] ~ Uniform(d.onset_lo_day[i] - 80.0, T_onset[i] - 1e-6)
            inc_i = T_onset[i] - T_inf[i]
            Turing.@addlogprob! logpdf(inc_dist, inc_i)
        else
            # Sourced case: T_inf anchored to listed exposure window.
            # GI > 0 enforced by rejecting trajectories where the secondary
            # was infected before its source.
            src = d.source_idx[i]
            T_inf[i] ~ Uniform(d.exp_lo_day[i], d.exp_hi_day[i])
            if T_inf[i] <= T_inf[src]
                Turing.@addlogprob! -Inf
            else
                inc_i  = T_onset[i] - T_inf[i]
                δ_pair = T_inf[i] - T_onset[src]
                Turing.@addlogprob! logpdf(inc_dist, inc_i)
                Turing.@addlogprob! logpdf(Normal(μ_δ, σ_δ), δ_pair)
            end
        end
        # Offspring count for case i (observed attributed secondaries).
        # log_R is clamped to a wide but finite range: under the non-centred
        # RW, `InitFromUniform` can produce an early-iterate log_R large
        # enough that `R = exp(log_R)` overflows and breaks the NB
        # `0 < p ≤ 1` check. The clamp is invisible once warmed up
        # (the posterior over log_R lives well inside [-50, 50]).
        R_i = exp(clamp(log_R[which_bin(T_inf[i], edges)], -50.0, 50.0))
        d.Zobs[i] ~ NegativeBinomial(k, k / (k + R_i))
    end
end
