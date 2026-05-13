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

@model function joint_model_def(d, edges)
    # Population-level parameters
    μ_inc ~ Normal(3.0, 0.5)                       # log-mean Inc (≈ log 20 d)
    σ_inc ~ truncated(Normal(0.0, 0.5); lower = 0) # log-SD Inc
    μ_δ   ~ Normal(0.0, 5.0)                       # population mean transmission timing (d from source onset)
    σ_δ   ~ truncated(Normal(0.0, 1.0); lower = 0) # population SD of transmission timing (d)
    # NB offspring dispersion via Stan's reciprocal-sqrt reparameterisation:
    # 1/√k is the SD multiplier in Var = μ + μ²·(1/√k)². Half-Normal(0, 1)
    # spans Poisson (1/√k → 0) to heavy super-spreader (1/√k ≈ 2)
    # symmetrically on the overdispersion scale.
    phi_inv_sqrt ~ truncated(Normal(0.0, 1.0); lower = 0)
    k := 1.0 / phi_inv_sqrt^2
    σ_rw  ~ truncated(Normal(0.0, 0.5); lower = 0) # log-R RW innovation SD allows sharp R(t) swings under interventions

    # Concrete element type derived from a sampled scalar — avoids the
    # dynamic-dispatch tax that `Vector{Real}` imposes on AD backends.
    T = typeof(μ_inc)

    # Non-centred random walk on log R(t) at the weekly knots. log_R[b] is
    # the value at knot b; R(t) is linearly interpolated between knots.
    n_knots = length(edges)
    log_R_init ~ Normal(log(1.5), 1.0)
    ε ~ Turing.filldist(Normal(zero(T), one(T)), n_knots - 1)
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
        # Clamp log R(t) to keep p = k/(k+R) strictly in (0, 1) during NUTS
        # exploration; the bounds sit well outside posterior support.
        R_i = exp(clamp(log_R_at(T_inf[i], edges, log_R), -50.0, 50.0))
        d.Zobs[i] ~ NegativeBinomial(k, k / (k + R_i))
    end
end
