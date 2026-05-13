## Joint model for the Epuyén ANDV outbreak.
##
## Three population-level components, written as swappable Turing models:
##   - Incubation period            — `incubation_model(μ_prior, σ_prior)`
##   - Transmission timing δ        — `transmission_delta_model(...)`
##   - log R(t) time series         — `random_walk_rt_model(n_bins; ...)`
##
## Per-case latents: an infection time `T_inf[i]` and an onset time
## `T_onset[i]`. The GI > 0 constraint `T_inf[secondary] > T_inf[source]`
## is enforced by a `-Inf` reject in the per-case loop.
##
## Real-time mode (gated on `d.obs_time !== nothing`):
##   - Index cases get Inc-only right-truncation `-logcdf(inc.dist, ...)`
##     (no δ).
##   - Sourced cases get the joint right-truncation
##     `-log F_offspring(obs_time − T_onset[src])`: the observation event
##     `T_onset[i] ≤ obs_time` is `δ + Inc(sec) ≤ obs_time − T_onset[src]`,
##     a joint constraint on the pair `(δ, Inc)` normalised by a single
##     `F_offspring`, not the product of marginal CDFs.
##   - NB rate `R` is thinned by the same `F_offspring(obs_time − T_onset[src])`
##     (binomial thinning of the offspring count). Same value as the
##     truncation, so a single vectorised `F_offspring` call serves both.

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
@model function incubation_model(μ_prior = Normal(3.0, 0.5),
                                 σ_prior = truncated(Normal(0.0, 0.5); lower = 0))
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
@model function transmission_delta_model(μ_prior = Normal(0.0, 5.0),
                                         σ_prior = truncated(Normal(0.0, 1.0); lower = 0))
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
    safe_nb(k, R)

`NegativeBinomial(k, p)` with `p = max(k/(k+R), eps(typeof(k)))`. The
clamp keeps the gradient finite when an extreme NUTS proposal overflows
`exp(log_R)` to `Inf`; the clamped value gives a vanishing likelihood
for any observed `Zobs > 0`, so the proposal is still rejected on
acceptance.
"""
safe_nb(k, R) = NegativeBinomial(k, max(k / (k + R), eps(typeof(k))))

"""
    joint_model(d, edges, foffspring_alg = _F_OFFSPRING_ALG;
                incubation   = incubation_model(),
                transmission = transmission_delta_model(),
                rt           = random_walk_rt_model(length(edges) + 1),
                k_prior      = truncated(Normal(0.3, 0.5); lower = 0))

Joint Bayesian model over incubation, transmission timing, and the
weekly random walk on log R(t). The three population components are
passed in as Turing submodels so priors and structural choices can be
swapped without editing this function — see [`incubation_model`](@ref),
[`transmission_delta_model`](@ref), and [`random_walk_rt_model`](@ref)
for the default contracts.
"""
@model function joint_model(d, edges, foffspring_alg = _F_OFFSPRING_ALG;
                            incubation   = incubation_model(),
                            transmission = transmission_delta_model(),
                            rt           = random_walk_rt_model(length(edges) + 1),
                            k_prior      = truncated(Normal(0.3, 0.5); lower = 0))
    inc    ~ to_submodel(incubation, false)
    delta  ~ to_submodel(transmission, false)
    _log_R ~ to_submodel(rt, false)
    k      ~ k_prior

    log_R := _log_R

    # Concrete element type so per-case latent vectors don't fall back
    # to `Vector{Real}` and tax AD with dynamic dispatch.
    T = typeof(inc.μ)

    T_onset = Vector{T}(undef, d.N)
    for i in 1:d.N
        T_onset[i] ~ Uniform(d.onset_lo_day[i], d.onset_hi_day[i])
    end

    realtime = d.obs_time !== nothing

    # F_offspring(obs_time − T_onset[src]) is the joint right-truncation
    # for an observed sourced pair and the binomial thinning factor for
    # source `src`'s offspring count — same probability, used in both
    # the per-case loop and the NB likelihood below.
    thins = realtime ?
        F_offspring(d.obs_time .- T_onset, inc.dist, delta.dist;
                    alg = foffspring_alg) : T[]

    T_inf = Vector{T}(undef, d.N)
    for i in 1:d.N
        src = d.source_idx[i]
        if src == 0
            T_inf[i] ~ Uniform(d.onset_lo_day[i] - 80.0, T_onset[i] - 1e-6)
            Turing.@addlogprob! logpdf(inc.dist, T_onset[i] - T_inf[i])
            realtime && Turing.@addlogprob! -logcdf(inc.dist, d.obs_time[i] - T_inf[i])
        else
            T_inf[i] ~ Uniform(d.exp_lo_day[i], d.exp_hi_day[i])
            if T_inf[i] <= T_inf[src]
                Turing.@addlogprob! oftype(zero(T), -Inf)
            else
                Turing.@addlogprob! logpdf(inc.dist,   T_onset[i] - T_inf[i])
                Turing.@addlogprob! logpdf(delta.dist, T_inf[i]   - T_onset[src])
                # `floatmin` guards Mooncake against a stray `log(0)`
                # when an extreme NUTS proposal collapses `thins[src]`.
                realtime && Turing.@addlogprob! -log(max(thins[src], floatmin(T)))
            end
        end
    end

    # `safe_nb` clamps the NB success probability into [eps, 1] so an
    # extreme NUTS proposal that overflows `exp(log_R)` to `Inf` doesn't
    # throw a DomainError on the gradient path.
    for i in 1:d.N
        R_i   = exp(log_R[which_bin(T_inf[i], edges)])
        R_eff = realtime ? R_i * thins[i] : R_i
        d.Zobs[i] ~ safe_nb(k, R_eff)
    end
end
