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
##   - Inc and δ contributions are right-truncated at the per-case cut-off
##     via `-logcdf` terms.
##   - The NB mean for source `src` is thinned by `F_cluster(obs_time[src] -
##     T_inf[src])` to account for offspring whose total chain has not yet
##     completed by the cut-off. When `obs_time === nothing` the whole
##     correction collapses out and the model is identical to the
##     retrospective form.
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
downstream quantities such as `F_cluster`.

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
    joint_model(d, edges, fcluster_alg = _F_CLUSTER_ALG;
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
distributions, the GI > 0 reject and the cluster-completeness thinning.
"""
@model function joint_model(d, edges, fcluster_alg = _F_CLUSTER_ALG;
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

    # Pass 1: sample T_inf and accrue the incubation / δ logprobs.
    # F_cluster is deferred to a single vectorised call so the 2-D
    # Gauss-Hermite tensor runs once across all cases instead of N times.
    T_inf = Vector{T}(undef, d.N)
    for i in 1:d.N
        if d.source_idx[i] == 0
            # Zoonotic index: free latent T_inf pre-onset.
            T_inf[i] ~ Uniform(d.onset_lo_day[i] - 80.0, T_onset[i] - 1e-6)
            inc_i = T_onset[i] - T_inf[i]
            Turing.@addlogprob! logpdf(inc.dist, inc_i)
            if realtime
                Turing.@addlogprob! -logcdf(inc.dist, d.obs_time[i] - T_inf[i])
            end
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
                Turing.@addlogprob! logpdf(inc.dist, inc_i)
                Turing.@addlogprob! logpdf(delta.dist, δ_pair)
                if realtime
                    Δ_inc = d.obs_time[i]   - T_inf[i]
                    Δ_δ   = d.obs_time[i]   - T_onset[src]
                    Turing.@addlogprob! -logcdf(inc.dist, Δ_inc)
                    Turing.@addlogprob! -logcdf(delta.dist, Δ_δ)
                end
            end
        end
    end

    # Pass 2: NB offspring likelihood. In realtime mode all N
    # cluster-completeness probabilities share the same population
    # distributions, so a single vector-valued F_cluster amortises the
    # quadrature across cases. Profiling showed F_cluster dominates the
    # per-eval cost; this collapses N quadrature solves into one.
    if realtime
        Δ_srcs = d.obs_time .- T_inf
        thins  = F_cluster(Δ_srcs, inc.dist, delta.dist; alg = fcluster_alg)
        for i in 1:d.N
            R_i   = exp(log_R[which_bin(T_inf[i], edges)])
            R_eff = R_i * thins[i]
            d.Zobs[i] ~ NegativeBinomial(k, k / (k + R_eff))
        end
    else
        for i in 1:d.N
            R_i = exp(log_R[which_bin(T_inf[i], edges)])
            d.Zobs[i] ~ NegativeBinomial(k, k / (k + R_i))
        end
    end
end
