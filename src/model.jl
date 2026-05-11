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
##   3. Time-varying reproduction   — log R(t) on a monthly random walk, with
##      number                        Negative-Binomial offspring (dispersion k).
##
## Each case has continuous latents: an infection time T_inf and a within-day
## onset offset. Daily-resolution onsets and exposure dates are handled by
## Bayesian data augmentation over these latents.
##
## The generation interval and serial interval are derived in post-processing
## from δ and Inc. GI ≥ 0 by definition; the per-pair constraint
## T_inf[secondary] > T_inf[source] enforces this at the latent level.

using Distributions, Turing

@model function joint_model(d, edges)
    # Population-level parameters
    μ_inc ~ Normal(3.0, 0.5)                                       # log-mean Inc (≈ log 20 d)
    σ_inc ~ truncated(Normal(0.0, 0.5); lower = 0.0, upper = 2.0)  # log-SD Inc
    μ_δ   ~ Normal(0.0, 5.0)                                       # population mean transmission timing (d from source onset)
    σ_δ   ~ truncated(Normal(0.0, 1.0); lower = 0.0, upper = 10.0) # population SD of transmission timing (d)
    k     ~ truncated(Normal(0.3, 0.5); lower = 0.0, upper = 10.0) # NB offspring dispersion (centred low — known super-spreader pathogen)
    σ_rw  ~ truncated(Normal(0.0, 0.5); lower = 0.0, upper = 2.0)  # log-R RW innovation SD

    # Random walk on log R(t) across the time bins
    n_bins = length(edges) + 1
    log_R = Vector{Real}(undef, n_bins)
    log_R[1] ~ Normal(log(1.5), 1.0)
    for b in 2:n_bins
        log_R[b] ~ Normal(log_R[b - 1], σ_rw)
    end

    inc_dist = LogNormal(μ_inc, σ_inc)

    # Latents per case: within-day onset offset and continuous infection time
    T_inf  = Vector{Real}(undef, d.N)
    offset = Vector{Real}(undef, d.N)
    for i in 1:d.N
        offset[i] ~ Uniform(0.0, 1.0)
    end
    T_onset = d.onset_day .+ offset

    for i in 1:d.N
        if d.source_idx[i] == 0
            # Zoonotic index: free latent T_inf pre-onset.
            T_inf[i] ~ Uniform(d.onset_day[i] - 80.0, T_onset[i] - 1e-6)
            inc_i = T_onset[i] - T_inf[i]
            Turing.@addlogprob! logpdf(inc_dist, inc_i)
        else
            # Sourced case: T_inf anchored to listed exposure window;
            # δ_pair contributes to (μ_δ, σ_δ). The hard constraint
            # T_inf[secondary] > T_inf[source] enforces GI > 0 per pair.
            src    = d.source_idx[i]
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
        # Offspring (per-case Z observed)
        R_i = exp(log_R[which_bin(T_inf[i], edges)])
        d.Zobs[i] ~ NegativeBinomial(k, k / (k + R_i))
    end
end
