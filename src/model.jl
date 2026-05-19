## Joint model for the Epuyén ANDV outbreak.
##
## Population-level components, written as swappable Turing submodels:
##   - Incubation period            — `incubation_model(...)`
##   - Transmission timing δ        — `transmission_delta_model(...)`
##   - log R(t) time series         — `random_walk_rt_model(n_knots; ...)`
##   - NB dispersion k              — `nb_dispersion_model(...)`
##
## Call chain:
##   joint_model
##     ├── delays_only_model
##     │     ├── incubation_model
##     │     ├── transmission_delta_model
##     │     └── latent_times_model
##     │           └── truncation_model   (real-time only)
##     └── case_model
##           ├── random_walk_rt_model
##           └── nb_dispersion_model

"""
$(TYPEDSIGNATURES)

Model sampling the location `μ_inc` and scale `σ_inc` of the incubation
period distribution. Returns `(; dist = dist_constructor(μ, σ), μ, σ)`.

The `dist_constructor` keyword lets the family be swapped without
editing the model. Defaults to `LogNormal`, in which case `μ`/`σ` are
the log-mean and log-SD of the incubation period.

# Keyword Arguments
- `dist_constructor`: two-argument distribution constructor called as
  `dist_constructor(μ, σ)`. Defaults to `LogNormal`.
- `μ_prior`: prior on the location parameter. Defaults to
  `Normal(3.0, 0.5)`.
- `σ_prior`: prior on the scale parameter, constrained positive.
  Defaults to `truncated(Normal(0.0, 0.5); lower = 0)`.
"""
@model function incubation_model(; dist_constructor = LogNormal,
        μ_prior = Normal(3.0, 0.5),
        σ_prior = truncated(Normal(0.0, 0.5); lower = 0))
    μ_inc ~ μ_prior
    σ_inc ~ σ_prior
    return (; dist = dist_constructor(μ_inc, σ_inc),
        μ = μ_inc, σ = σ_inc)
end

"""
$(TYPEDSIGNATURES)

Model for the population mean `μ_δ` and SD `σ_δ` of the per-pair
transmission timing. Returns `(; dist = dist_constructor(μ, σ), μ, σ)`.

# Keyword Arguments
- `dist_constructor`: two-argument distribution constructor called as
  `dist_constructor(μ, σ)`. Defaults to `Normal`.
- `μ_prior`: prior on the population mean `μ_δ`. Defaults to
  `Normal(0.0, 5.0)`.
- `σ_prior`: prior on the population SD `σ_δ`, constrained positive.
  Defaults to `truncated(Normal(0.0, 1.0); lower = 0)`.
"""
@model function transmission_delta_model(; dist_constructor = Normal,
        μ_prior = Normal(0.0, 5.0),
        σ_prior = truncated(Normal(0.0, 1.0); lower = 0))
    μ_δ ~ μ_prior
    σ_δ ~ σ_prior
    return (; dist = dist_constructor(μ_δ, σ_δ), μ = μ_δ, σ = σ_δ)
end

"""
$(TYPEDSIGNATURES)

Non-centred weekly random walk on log R(t) at `n_knots` knots.
Returns the length-`n_knots` `log_R` vector evaluated at the knot dates;
[`log_R_at`](@ref) linearly interpolates between knots.

# Arguments
- `n_knots`: number of weekly knot points at which `log_R` is evaluated.

# Keyword Arguments
- `init_prior`: prior on the initial log R(t) value `log_R_init`. Defaults
  to `Normal(log(1.5), 1.0)`.
- `sigma_prior`: prior on the random walk step SD `σ_rw`, constrained
  positive. Defaults to `truncated(Normal(0.0, 0.2); lower = 0)`.
"""
@model function random_walk_rt_model(n_knots::Integer;
        init_prior = Normal(log(1.5), 1.0),
        sigma_prior = truncated(Normal(0.0, 0.2); lower = 0))
    σ_rw ~ sigma_prior
    log_R_init ~ init_prior
    T = typeof(log_R_init)
    ε ~ Turing.filldist(Normal(zero(T), one(T)), n_knots - 1)
    log_R := vcat(log_R_init, log_R_init .+ accumulate(+, σ_rw .* ε))
    return (; log_R)
end

"""
$(TYPEDSIGNATURES)

Model for the negative binomial dispersion `k` using the Stan-default
`1/√k` reparameterisation. Samples `phi_inv_sqrt` from `phi_prior` and
returns `(; k = 1 / phi_inv_sqrt^2, phi_inv_sqrt)`.

# Arguments
- `phi_prior`: prior on `1/√k`, constrained positive. Defaults to
  `truncated(Normal(0.0, 1.0); lower = 0)`.
"""
@model function nb_dispersion_model(
        phi_prior = truncated(Normal(0.0, 1.0); lower = 0))
    phi_inv_sqrt ~ phi_prior
    k := 1.0 / phi_inv_sqrt^2
    return (; k, phi_inv_sqrt)
end

"""
    safe_nb(k, R)

`NegativeBinomial(k, p)` with `p = max(k/(k+R), eps(typeof(k)))`. The
clamp stabilises the gradient when an extreme NUTS proposal makes
`R` so large that `p` underflows to zero, which would otherwise leave
`logpdf` at `-Inf` and break AD.
"""
safe_nb(k, R) = NegativeBinomial(k, max(k / (k + R), eps(typeof(k))))

"""
$(TYPEDSIGNATURES)

Per-case negative binomial likelihood for offspring counts `Z`. Draws
`log_R` from the nested `rt` submodel and the NB dispersion `k` from
the nested `dispersion` submodel; both are surfaced as deterministic
chain variables via `:=`. In real-time mode, the per-case rate is
thinned by the offspring-completeness `p` (pass empty `p` for
retrospective mode, where `R_eff = R`).

# Arguments
- `Z`: vector of observed offspring counts per case.
- `edges`: vector of knot dates at which `log_R` is defined.
- `T_onset`: vector of per-case onset times used to index `log_R`.
- `p`: vector of per-case offspring-completeness values, or empty for
  retrospective mode.

# Keyword Arguments
- `rt`: Turing submodel for the log R(t) random walk. Defaults to
  `random_walk_rt_model(length(edges))`.
- `dispersion`: Turing submodel for the NB dispersion `k`. Defaults to
  `nb_dispersion_model()`.
"""
@model function case_model(Z, edges, T_onset, p;
        rt = random_walk_rt_model(length(edges)),
        dispersion = nb_dispersion_model())
    random_walk ~ to_submodel(rt, false)
    nb_dispersion ~ to_submodel(dispersion, false)
    realtime = !isempty(p)
    for i in eachindex(Z)
        R_i = exp(log_R_at(T_onset[i], edges, random_walk.log_R))
        R_eff = realtime ? R_i * p[i] : R_i
        Z[i] ~ safe_nb(nb_dispersion.k, R_eff)
    end
end

"""
$(TYPEDSIGNATURES)

Posterior-predictive draw for a single case under the [`case_model`](@ref)
negative-binomial likelihood. Uses the Gamma-Poisson conjugate structure
of NB: `λ ~ Gamma(k + Z_obs, R / (k + R · p_i))`, then
`Poisson(λ · prob)`. Real-time predictors dispatch on the joint
model's observation-submodel constructor (`model.defaults.cases`) so
alternative obs models (e.g. Poisson, beta-binomial) can supply their
own `posterior_predictive` methods without editing the predictor.

# Arguments
- `cases_sub`: the observation-submodel constructor stored on
  `joint_model`'s defaults; here `typeof(case_model)`.
- `rng`: RNG used for the Gamma and Poisson draws.
- `k`: per-draw NB dispersion `k`.
- `Z_obs`: observed offspring count for this source.
- `R`: per-source rate `R_i = exp(log_R_at(T_onset_i, edges, log_R))`.
- `p_i`: offspring-completeness probability for this source.
- `prob`: per-source future-onset thinning probability (e.g.
  `q_i − p_i` for the controlled counterfactual).
"""
function posterior_predictive(::typeof(case_model),
        rng, k, Z_obs, R, p_i, prob)
    λ = rand(rng, Gamma(k + Z_obs, R / (k + R * p_i)))
    return rand(rng, Poisson(λ * prob))
end

"""
$(TYPEDSIGNATURES)

Real-time right-truncation submodel. Adds two contributions to the
log-likelihood:

1. Per-case Inc right-truncation `-logcdf(inc_dist, obs_time − T_inf[i])`
   for **index cases only**. Sourced cases already get their offspring-
   chain right-truncation through the convolved-delay denominator in
   contribution (2), so applying `-logcdf` to them would double-count.
2. Per-pair offspring-completeness denominator `-log(p[src])` for each
   sourced case, with
   `p = cdf(ConvolvedDelays(inc_dist, delta_dist), obs_time .- T_onset)`.

Both contributions are vectorised into a single `@addlogprob!` call
each rather than per-case loops. Returns `(; p)`, where `p` is the
offspring-chain completion cdf
`cdf(ConvolvedDelays(inc_dist, delta_dist), obs_time .- T_onset)`.
Called with empty `T_inf` / `T_onset` / `source_idx` in retrospective
mode so both `@addlogprob!` calls reduce to `0` and `p` comes back
empty.

# Arguments
- `T_inf`: per-case infection times.
- `T_onset`: per-case onset times.
- `source_idx`: per-case source indices (0 for index cases).
- `obs_time`: per-case real-time cut-off as a vector of day numbers,
  same length as `T_inf` and `T_onset`.
- `inc_dist`: incubation period distribution.
- `delta_dist`: per-pair transmission timing distribution.
"""
@model function truncation_model(T_inf, T_onset, source_idx, obs_time,
        inc_dist, delta_dist)
    T = eltype(T_onset)
    index = findall(==(0), source_idx)
    Turing.@addlogprob! -sum(logcdf.(inc_dist,
        obs_time[index] .- T_inf[index]))
    p = cdf(ConvolvedDelays(inc_dist, delta_dist), obs_time .- T_onset)
    sourced = source_idx[source_idx .!= 0]
    Turing.@addlogprob! -sum(log.(max.(p[sourced], floatmin(T))))
    return (; p)
end

"""
$(TYPEDSIGNATURES)

Per-case latent infection and onset time submodel. Samples `T_onset[i]`
and `T_inf[i]` for each case and adds the marginal incubation and
transmission-timing log-densities for the resulting pairs. Sourced
cases enforce the GI > 0 constraint by rejecting infeasible draws.

The real-time right-truncation contributions live in the nested
[`truncation_model`](@ref), invoked unconditionally with empty inputs
in retrospective mode so the truncation reduces to a no-op.

# Arguments
- `d`: structured line-list data as returned by `build_data`, providing
  per-case onset bounds, exposure bounds, source indices, and optional
  `obs_time` for real-time fits.
- `inc_dist`: incubation period distribution from the incubation submodel.
- `δ_dist`: transmission-timing distribution from the transmission submodel.
"""
@model function latent_times_model(d, inc_dist, δ_dist)
    T = partype(inc_dist)
    T_onset = Vector{T}(undef, d.N)
    for i in 1:d.N
        T_onset[i] ~ Uniform(d.onset_lo_day[i], d.onset_hi_day[i])
    end
    T_inf = Vector{T}(undef, d.N)
    for i in 1:d.N
        src = d.source_idx[i]
        if src == 0
            T_inf[i] ~ Uniform(d.onset_lo_day[i] - 80.0,
                T_onset[i] - 1e-6)
            Turing.@addlogprob! logpdf(inc_dist, T_onset[i] - T_inf[i])
        else
            T_inf[i] ~ Uniform(d.exp_lo_day[i],
                min(d.exp_hi_day[i], T_onset[i] - 1e-6))
            if T_inf[i] <= T_inf[src]
                Turing.@addlogprob! oftype(zero(T), -Inf)
            else
                Turing.@addlogprob! logpdf(inc_dist,
                    T_onset[i] - T_inf[i])
                Turing.@addlogprob! logpdf(δ_dist,
                    T_inf[i] - T_onset[src])
            end
        end
    end

    if d.obs_time !== nothing
        truncation ~ to_submodel(
            truncation_model(T_inf, T_onset, d.source_idx,
                d.obs_time, inc_dist, δ_dist), false)
        p = truncation.p
    else
        p = T[]
    end

    return (; T_inf, T_onset, p)
end

"""
$(TYPEDSIGNATURES)

Delays-only model: composes the incubation, transmission timing, and
per-case latent submodels into a standalone diagnostic fit. Useful
when the full joint fit collapses: a clean delays-only fit isolates
any pathology to the R(t) / `case_model` half of the joint likelihood
rather than the delay parameters.

Returns `(; T_inf, T_onset, p, inc, δ)`. The `p` field carries the
per-case offspring-completeness vector from
[`latent_times_model`](@ref) (empty in retrospective mode); `inc` and
`δ` are the realised incubation and transmission timing distributions
so downstream submodels can re-use them.

# Arguments
- `d`: structured line-list data as returned by `build_data`,
  including per-case onset bounds, exposure bounds, source indices,
  and optional `obs_time` for real-time fits.

# Keyword Arguments
- `incubation`: Turing submodel for the incubation period. Defaults to
  `incubation_model()`.
- `transmission`: Turing submodel for the transmission timing `δ`.
  Defaults to `transmission_delta_model()`.
"""
@model function delays_only_model(d;
        incubation = incubation_model(),
        transmission = transmission_delta_model())
    inc ~ to_submodel(incubation, false)
    delta ~ to_submodel(transmission, false)
    latent ~ to_submodel(
        latent_times_model(d, inc.dist, delta.dist), false)
    return (; T_inf = latent.T_inf, T_onset = latent.T_onset,
        p = latent.p, inc = inc.dist, δ = delta.dist)
end

"""
$(TYPEDSIGNATURES)

Joint Bayesian model. Composes [`delays_only_model`](@ref) (incubation
+ δ + per-case latents + real-time right-truncation) with
[`case_model`](@ref) (R(t) random walk + NB dispersion + per-case NB
offspring likelihood). The case-level outputs `log_R` and `k` are
surfaced via `:=` from inside `case_model`.

# Arguments
- `d`: structured line-list data as returned by `build_data`, including
  per-case onset bounds, exposure bounds, source indices, offspring counts
  `Zobs`, and optional `obs_time` for real-time fits.
- `edges`: vector of knot dates (as day numbers) at which `log_R` is
  defined.

# Keyword Arguments
- `incubation`: Turing submodel for the incubation period passed
  through to `delays_only_model`. Defaults to `incubation_model()`.
- `transmission`: Turing submodel for the transmission timing `δ`
  passed through to `delays_only_model`. Defaults to
  `transmission_delta_model()`.
- `rt`: Turing submodel for the log R(t) random walk passed through to
  `case_model`. Defaults to `random_walk_rt_model(length(edges))`.
- `dispersion`: Turing submodel for the NB dispersion `k` passed
  through to `case_model`. Defaults to `nb_dispersion_model()`.
- `cases`: observation-submodel constructor invoked as
  `cases(d.Zobs, edges, delays.T_onset, delays.p; rt, dispersion)`.
  Defaults to [`case_model`](@ref). Exposed on `model.defaults.cases`
  so real-time predictors can dispatch [`posterior_predictive`](@ref)
  on the obs submodel.
"""
@model function joint_model(d, edges;
        incubation = incubation_model(),
        transmission = transmission_delta_model(),
        rt = random_walk_rt_model(length(edges)),
        dispersion = nb_dispersion_model(),
        cases = case_model)
    delays ~ to_submodel(
        delays_only_model(d; incubation, transmission), false)
    case_obs ~ to_submodel(
        cases(d.Zobs, edges, delays.T_onset, delays.p;
            rt, dispersion), false)
end
