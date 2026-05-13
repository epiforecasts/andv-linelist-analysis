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
##
## Real-time mode (gated on `d.obs_time !== nothing`):
##   - Index (zoonotic) cases: Inc is right-truncated at the per-case
##     cut-off via a single `-logcdf(inc.dist, obs_time - T_inf[i])` term.
##   - Sourced cases: the observation event T_onset[i] ≤ obs_time is a
##     joint constraint on the offspring pair `(δ, Inc)` —
##     `δ + Inc ≤ obs_time − T_onset[src]` — and so is normalised by a
##     single `-log F_offspring(obs_time − T_onset[src])`, not the product
##     of marginal Inc and δ CDFs.
##   - The NB mean for source `src` is thinned by `F_offspring(obs_time[src] -
##     T_onset[src])` to account for offspring whose `δ + Inc(sec)` chain
##     has not yet completed by the cut-off. The source's own incubation is
##     already pinned by the sampled latents `T_onset[src]` and `T_inf[src]`,
##     so only `δ` and `Inc(sec)` are marginalised here. When
##     `obs_time === nothing` the whole correction collapses out and the
##     model is identical to the retrospective form.
##
## Structure: the three population-level components are written as Turing
## models, each parameterised by the priors it samples from. Swapping a
## prior — or, by writing a parallel model with the same return contract,
## swapping the distributional family or the R(t) time-series model — does
## not require touching `joint_model`.

"""
    incubation_model(μ_prior, σ_prior)

Model sampling the log-mean `μ_inc` and log-SD `σ_inc` of a LogNormal
incubation period. Returns a `NamedTuple` `(dist, μ, σ)` where `dist`
is `LogNormal(μ_inc, σ_inc)`. The parent model uses `dist` to score
per-case incubation contributions and the raw parameters to evaluate
downstream quantities such as `F_offspring`.

The priors are arguments so they can be changed without editing this
file. To swap the *family* (Gamma, Weibull, …) write a parallel
model that returns a NamedTuple with the same fields.
"""
@model function incubation_model(μ_prior, σ_prior)
    μ_inc ~ μ_prior
    σ_inc ~ σ_prior
    return (; dist = LogNormal(μ_inc, σ_inc), μ = μ_inc, σ = σ_inc)
end

"""
    transmission_delta_model(μ_prior, σ_prior)

Model for the population mean `μ_δ` and SD `σ_δ` of the per-pair
transmission timing (gap between a secondary's infection and its
source's onset). Returns `(dist = Normal(μ_δ, σ_δ), μ, σ)`. As with
[`incubation_model`](@ref) the priors are arguments and the return
contract is the swap point for alternative families.
"""
@model function transmission_delta_model(μ_prior, σ_prior)
    μ_δ ~ μ_prior
    σ_δ ~ σ_prior
    return (; dist = Normal(μ_δ, σ_δ), μ = μ_δ, σ = σ_δ)
end

"""
    random_walk_rt_model(n_bins;
                            init_prior  = Normal(log(1.5), 1.0),
                            sigma_prior = truncated(Normal(0.0, 0.5); lower = 0))

Non-centred weekly random walk on log R(t) over `n_bins` bins. Returns
the length-`n_bins` `log_R` vector. The non-centred parameterisation
decouples `σ_rw` from the walk to avoid the funnel that diverges NUTS
under the centred form.

`joint_model` accepts any model that returns a length-`n_bins`
real-valued vector for log R(t), so alternative time-series structures
(AR1, GP, piecewise constant) drop in by writing a parallel model with
the same return contract.
"""
@model function random_walk_rt_model(n_bins::Integer;
                                        init_prior  = Normal(log(1.5), 1.0),
                                        sigma_prior = truncated(Normal(0.0, 0.5); lower = 0))
    σ_rw       ~ sigma_prior
    log_R_init ~ init_prior
    T = typeof(log_R_init)
    ε ~ Turing.filldist(Normal(zero(T), one(T)), n_bins - 1)
    return vcat(log_R_init, log_R_init .+ accumulate(+, σ_rw .* ε))
end

"""
    joint_model(d, edges, foffspring_alg = _F_OFFSPRING_ALG;
                incubation   = incubation_model(Normal(3.0, 0.5),
                                                   truncated(Normal(0.0, 0.5); lower = 0)),
                transmission = transmission_delta_model(Normal(0.0, 5.0),
                                                           truncated(Normal(0.0, 1.0); lower = 0)),
                rt           = random_walk_rt_model(length(edges) + 1),
                k_prior      = truncated(Normal(0.3, 0.5); lower = 0))

Joint Bayesian model over incubation, transmission timing and the
weekly random walk on log R(t). The three population-level components
are passed in as Turing submodels so priors and structural choices can
be swapped without editing this function. See
[`incubation_model`](@ref), [`transmission_delta_model`](@ref)
and [`random_walk_rt_model`](@ref) for the default contracts.

The per-case incubation / δ / NB contributions are kept inline because
factoring them into submodels would shuffle complexity rather than
remove it: the two source-vs-index branches share the population
distributions, the GI > 0 reject and the offspring-completeness thinning.
"""
@model function joint_model(d, edges, foffspring_alg = _F_OFFSPRING_ALG;
                            incubation   = incubation_model(
                                Normal(3.0, 0.5),
                                truncated(Normal(0.0, 0.5); lower = 0)),
                            transmission = transmission_delta_model(
                                Normal(0.0, 5.0),
                                truncated(Normal(0.0, 1.0); lower = 0)),
                            rt           = random_walk_rt_model(length(edges) + 1),
                            k_prior      = truncated(Normal(0.3, 0.5); lower = 0))
    inc    ~ to_submodel(incubation, false)
    delta  ~ to_submodel(transmission, false)
    _log_R ~ to_submodel(rt, false)
    k      ~ k_prior

    # Track log_R as a deterministic so MCMCChains sees the same `log_R`
    # vector the pre-refactor model exposed.
    log_R := _log_R

    # Concrete element type derived from a sampled scalar — avoids the
    # dynamic-dispatch tax that `Vector{Real}` imposes on AD backends.
    T = typeof(inc.μ)

    # T_onset is a latent over the recorded onset window (defaults to a
    # one-day window when only a single onset date was recorded).
    T_onset = Vector{T}(undef, d.N)
    for i in 1:d.N
        T_onset[i] ~ Uniform(d.onset_lo_day[i], d.onset_hi_day[i])
    end

    realtime = d.obs_time !== nothing

    # In realtime mode the same `F_offspring(obs_time − T_onset[src])`
    # appears in two places — the NB thinning of source `src` (pass 2)
    # and the joint right-truncation of every offspring whose recorded
    # source is `src` (pass 1, sourced branch). Compute it once across
    # all N cases via a single vector-valued solve so the 1-D
    # Gauss-Hermite rule runs once instead of N+M times; both passes
    # then index into the resulting vector.
    thins = if realtime
        F_offspring(d.obs_time .- T_onset, inc.dist, delta.dist;
                    alg = foffspring_alg)
    else
        T[]
    end

    # Pass 1: sample T_inf and accrue the incubation / δ logprobs.
    # Per-case logprob contributions are accumulated locally into `lp`
    # and added in a single `@addlogprob!` at the loop end. Each
    # `@addlogprob!` goes through the DynamicPPL accumulator machinery
    # (and Mooncake records each as a separate tape entry); collapsing
    # the 1-to-4 addlogprob calls per case into one removes most of
    # that overhead.
    T_inf = Vector{T}(undef, d.N)
    lp = zero(T)
    for i in 1:d.N
        if d.source_idx[i] == 0
            # Zoonotic index: free latent T_inf pre-onset.
            T_inf[i] ~ Uniform(d.onset_lo_day[i] - 80.0, T_onset[i] - 1e-6)
            inc_i = T_onset[i] - T_inf[i]
            lp += logpdf(inc.dist, inc_i)
            if realtime
                lp -= logcdf(inc.dist, d.obs_time[i] - T_inf[i])
            end
        else
            # Sourced case: T_inf anchored to listed exposure window.
            # GI > 0 enforced by rejecting trajectories where the secondary
            # was infected before its source.
            src = d.source_idx[i]
            T_inf[i] ~ Uniform(d.exp_lo_day[i], d.exp_hi_day[i])
            if T_inf[i] <= T_inf[src]
                lp += oftype(lp, -Inf)
            else
                inc_i  = T_onset[i] - T_inf[i]
                δ_pair = T_inf[i] - T_onset[src]
                lp += logpdf(inc.dist, inc_i)
                lp += logpdf(delta.dist, δ_pair)
                if realtime
                    # Joint truncation on the offspring pair (δ, Inc):
                    # the observation event T_onset[i] ≤ obs_time is
                    # equivalent to δ + Inc(sec) ≤ obs_time − T_onset[src],
                    # i.e. a single F_offspring normaliser rather than the
                    # product of marginal Inc and δ CDFs. Reuses the
                    # pre-computed `thins[src]` so the quadrature is
                    # shared with the pass-2 NB thinning. `floatmin`
                    # guards against `log(0)` when `Δ_srcs[src]` lands
                    # at or below zero on a NUTS step (T_onset[src] at
                    # the upper end of its onset window) — the AD path
                    # would otherwise hit a DomainError under Mooncake.
                    lp -= log(max(thins[src], floatmin(T)))
                end
            end
        end
    end
    Turing.@addlogprob! lp

    # Pass 2: NB offspring likelihood. The reference point for thinning
    # is the source's onset time (a sampled latent), not its infection
    # time: the source's own incubation is already scored, so the
    # remaining offspring delay is `δ + Inc(sec)`. `thins` was computed
    # above so it could be re-used for the per-offspring truncation.
    #
    # `nb_p` clamps the NB success probability into [eps, 1] so that an
    # extreme NUTS proposal pushing `log_R[bin]` past ~700 (overflow to
    # Inf in R_i) does not throw a `DomainError` inside the
    # differentiated path. The clamped value still produces a vanishing
    # NB likelihood for any observed `Zobs > 0`, so the proposal is
    # rejected at acceptance — the clamp just keeps the gradient finite.
    nb_p(k, R) = max(k / (k + R), eps(typeof(k)))
    if realtime
        for i in 1:d.N
            R_i   = exp(log_R[which_bin(T_inf[i], edges)])
            R_eff = R_i * thins[i]
            d.Zobs[i] ~ NegativeBinomial(k, nb_p(k, R_eff))
        end
    else
        for i in 1:d.N
            R_i = exp(log_R[which_bin(T_inf[i], edges)])
            d.Zobs[i] ~ NegativeBinomial(k, nb_p(k, R_i))
        end
    end
end
