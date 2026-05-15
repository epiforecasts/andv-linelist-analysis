## Joint model for the Epuyén ANDV outbreak.
##
## Three population-level components, written as swappable Turing models:
##   - Incubation period            — `incubation_model(μ_prior, σ_prior)`
##   - Transmission timing δ        — `transmission_delta_model(...)`
##   - log R(t) time series         — `random_walk_rt_model(n_knots; ...)`
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
$(TYPEDSIGNATURES)

Model sampling the log-mean `μ_inc` and log-SD `σ_inc` of a LogNormal
incubation period.
Returns `(dist = LogNormal(μ_inc, σ_inc), μ, σ)`.

# Arguments
- `μ_prior`: prior on the log-mean `μ_inc`. Defaults to `Normal(3.0, 0.5)`.
- `σ_prior`: prior on the log-SD `σ_inc`, constrained positive. Defaults to
  `truncated(Normal(0.0, 0.5); lower = 0)`.
"""
@model function incubation_model(μ_prior = Normal(3.0, 0.5),
        σ_prior = truncated(Normal(0.0, 0.5); lower = 0))
    μ_inc ~ μ_prior
    σ_inc ~ σ_prior
    return (; dist = LogNormal(μ_inc, σ_inc), μ = μ_inc, σ = σ_inc)
end

"""
$(TYPEDSIGNATURES)

Model for the population mean `μ_δ` and SD `σ_δ` of the per-pair
transmission timing.
Returns `(dist = Normal(μ_δ, σ_δ), μ, σ)`.

# Arguments
- `μ_prior`: prior on the population mean `μ_δ`. Defaults to
  `Normal(0.0, 5.0)`.
- `σ_prior`: prior on the population SD `σ_δ`, constrained positive.
  Defaults to `truncated(Normal(1.0, 0.5); lower = 0)`.
"""
@model function transmission_delta_model(μ_prior = Normal(0.0, 5.0),
        σ_prior = truncated(Normal(1.0, 0.5); lower = 0))
    μ_δ ~ μ_prior
    σ_δ ~ σ_prior
    return (; dist = Normal(μ_δ, σ_δ), μ = μ_δ, σ = σ_δ)
end

"""
$(TYPEDSIGNATURES)

Non-centred weekly random walk on log R(t) at `n_knots` knots.
Returns the length-`n_knots` `log_R` vector evaluated at the knot dates;
`log_R_at` linearly interpolates between knots.

# Arguments
- `n_knots`: number of weekly knot points at which `log_R` is evaluated.

# Keyword Arguments
- `init_prior`: prior on the initial log R(t) value `log_R_init`. Defaults
  to `Normal(log(1.5), 0.5)`.
- `sigma_prior`: prior on the random walk step SD `σ_rw`, constrained
  positive. Defaults to `truncated(Normal(0.0, 0.5); lower = 0)`.
"""
@model function random_walk_rt_model(n_knots::Integer;
        init_prior = Normal(log(1.5), 0.5),
        sigma_prior = truncated(Normal(0.0, 0.5); lower = 0))
    σ_rw ~ sigma_prior
    log_R_init ~ init_prior
    T = typeof(log_R_init)
    ε ~ Turing.filldist(Normal(zero(T), one(T)), n_knots - 1)
    return vcat(log_R_init, log_R_init .+ accumulate(+, σ_rw .* ε))
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
    return (; k = 1.0 / phi_inv_sqrt^2, phi_inv_sqrt)
end

"""
    safe_nb(k, R)

`NegativeBinomial(k, p)` with `p = max(k/(k+R), eps(typeof(k)))`. Keeps
the gradient finite when an extreme NUTS proposal overflows `exp(log_R)`
to `Inf`.
"""
safe_nb(k, R) = NegativeBinomial(k, max(k / (k + R), eps(typeof(k))))

"""
$(TYPEDSIGNATURES)

Real-time right-truncation submodel. Evaluates the per-case
offspring-completeness `p_i = cdf(delay_dist, obs_time − T_onset[i])`
under the [`CombinedDelay`](@ref) `δ + Inc(sec)` and adds the per-pair
thinning denominator `-log(p[src])` once per sourced case, completing
the joint right-truncation correction for observed transmission pairs.
Returns `(; p)`.

The vector of `p` values is computed via a single
[`F_offspring`](@ref) quadrature solve rather than per-element
`cdf(delay_dist, ·)` calls, so AD passes are O(quadrature nodes)
rather than O(N × quadrature nodes).

# Arguments
- `T_onset`: per-case onset times.
- `source_idx`: per-case source indices (0 for index cases).
- `obs_time`: real-time cut-off as a day number.
- `delay_dist`: [`CombinedDelay`](@ref) of `δ + Inc(sec)`.

# Keyword Arguments
- `foffspring_alg`: quadrature algorithm passed to `F_offspring`.
  Defaults to `_F_OFFSPRING_ALG`.
"""
@model function truncation_model(T_onset, source_idx, obs_time,
        delay_dist::CombinedDelay;
        foffspring_alg = _F_OFFSPRING_ALG)
    T = eltype(T_onset)
    Δ = obs_time .- T_onset
    p = F_offspring(Δ, delay_dist.inc, delay_dist.δ; alg = foffspring_alg)
    for i in eachindex(source_idx)
        src = source_idx[i]
        src == 0 && continue
        Turing.@addlogprob! -log(max(p[src], floatmin(T)))
    end
    return (; p)
end

"""
$(TYPEDSIGNATURES)

Per-case negative binomial likelihood for offspring counts `Z`, thinned
by the per-case offspring-completeness `p` in real-time mode. Pass an
empty `p` for retrospective mode (`R_eff = R`).

# Arguments
- `Z`: vector of observed offspring counts per case.
- `T_onset`: vector of per-case onset times used to index `log_R`.
- `edges`: vector of knot dates at which `log_R` is defined.
- `log_R`: vector of log R(t) values at the knot dates.
- `k`: negative binomial dispersion.
- `p`: vector of per-case offspring-completeness values, or empty for
  retrospective mode.
"""
@model function case_model(Z, T_onset, edges, log_R, k, p)
    realtime = !isempty(p)
    for i in eachindex(Z)
        R_i = exp(log_R_at(T_onset[i], edges, log_R))
        R_eff = realtime ? R_i * p[i] : R_i
        Z[i] ~ safe_nb(k, R_eff)
    end
end

"""
$(TYPEDSIGNATURES)

Per-case latent infection and onset time submodel. Samples `T_onset[i]`
and `T_inf[i]` for each case and adds the marginal incubation and
transmission-timing log-densities for the resulting pairs. Index cases
get the Inc right-truncation under `obs_time` in real-time mode; sourced
cases enforce the GI > 0 constraint by rejecting infeasible draws.

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
    realtime = d.obs_time !== nothing
    T_inf = Vector{T}(undef, d.N)
    for i in 1:d.N
        src = d.source_idx[i]
        if src == 0
            T_inf[i] ~ Uniform(d.onset_lo_day[i] - 80.0, T_onset[i] - 1e-6)
            Turing.@addlogprob! logpdf(inc_dist, T_onset[i] - T_inf[i])
            realtime &&
                Turing.@addlogprob! -logcdf(inc_dist, d.obs_time[i] - T_inf[i])
        else
            T_inf[i] ~ Uniform(d.exp_lo_day[i],
                min(d.exp_hi_day[i], T_onset[i] - 1e-6))
            if T_inf[i] <= T_inf[src]
                Turing.@addlogprob! oftype(zero(T), -Inf)
            else
                Turing.@addlogprob! logpdf(inc_dist, T_onset[i] - T_inf[i])
                Turing.@addlogprob! logpdf(δ_dist, T_inf[i] - T_onset[src])
            end
        end
    end
    return (; T_inf, T_onset)
end

"""
$(TYPEDSIGNATURES)

Delays-only joint model: fits the incubation period and transmission
timing submodels against the augmented line-list latents `T_inf` and
`T_onset` but omits the R(t) random walk, the NB offspring likelihood,
and the offspring-completeness thinning correction.

Useful as both a standalone diagnostic fit and the first submodel of
[`joint_model`](@ref): a delays-only fit that converges where the full
joint fit collapses isolates the pathology to the R(t) / `case_model`
half of the joint likelihood rather than the delay parameters. Index
cases retain the same Inc right-truncation under `obs_time` as in the
full model; sourced cases get the marginal Inc and δ log-densities
with the GI > 0 constraint enforced, but no joint `F_offspring`
truncation (pure delay fitting).

Returns `(; T_inf, T_onset, inc, δ)` so the joint model can pass the
incubation and δ distributions through to the truncation and case
submodels without re-instantiating them.

# Arguments
- `d`: structured line-list data as returned by `build_data`,
  including per-case onset bounds, exposure bounds, source indices, and
  optional `obs_time` for real-time fits.

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
        inc = inc.dist, δ = delta.dist)
end

"""
$(TYPEDSIGNATURES)

Joint Bayesian model over incubation, transmission timing, and the
weekly random walk on log R(t) at the knot dates.
The delays bundle (incubation + δ + latent times), R(t), and
dispersion components are all passed in as Turing submodels so priors
and structural choices can be swapped without editing this function.

The submodel sequence is: `delays` → (real-time only) `truncation` →
`rt` + `dispersion` → `cases`. The truncation submodel is invoked
unconditionally and given empty inputs in retrospective mode so its
`p` vector comes back empty and `case_model` skips the thinning.

# Arguments
- `d`: structured line-list data as returned by `build_data`, including
  per-case onset bounds, exposure bounds, source indices, offspring counts
  `Zobs`, and optional `obs_time` for real-time fits.
- `edges`: vector of knot dates (as day numbers) at which `log_R` is
  defined.
- `foffspring_alg`: quadrature algorithm passed to `F_offspring` for the
  real-time offspring-completeness integral.

# Keyword Arguments
- `delays`: Turing submodel bundling incubation, transmission timing
  and the per-case latents. Defaults to `delays_only_model(d)`.
- `rt`: Turing submodel for the log R(t) random walk. Defaults to
  `random_walk_rt_model(length(edges))`.
- `dispersion`: Turing submodel for the NB dispersion `k` (via the
  Stan-default `1/√k` reparameterisation). Defaults to
  `nb_dispersion_model()`; swap to plug in a different `phi_prior`.
"""
@model function joint_model(d, edges, foffspring_alg = _F_OFFSPRING_ALG;
        delays = delays_only_model(d),
        rt = random_walk_rt_model(length(edges)),
        dispersion = nb_dispersion_model())
    dly ~ to_submodel(delays, false)
    T_onset = dly.T_onset
    delay_dist = CombinedDelay(dly.inc, dly.δ)

    realtime = d.obs_time !== nothing
    T = eltype(T_onset)
    trunc_T_onset = realtime ? T_onset : T[]
    trunc_src = realtime ? d.source_idx : Int[]
    obs_time_t = realtime ? d.obs_time : T[]
    trunc ~ to_submodel(
        truncation_model(trunc_T_onset, trunc_src, obs_time_t,
            delay_dist; foffspring_alg = foffspring_alg), false)

    _log_R ~ to_submodel(rt, false)
    disp ~ to_submodel(dispersion, false)
    log_R := _log_R
    k := disp.k

    cases ~ to_submodel(
        case_model(d.Zobs, T_onset, edges, log_R, k, trunc.p), false)
end
