## Real-time extension utilities.

import FastGaussQuadrature
using Integrals: IntegralProblem, GaussLegendre, solve

"""
$(TYPEDSIGNATURES)

Trim the line-list to cases whose onset is on or before `obs_date`
(those an analyst would know about at the cut-off).
Source attributions pointing outside the retained set are dropped and
`Z` is rebuilt to count only retained offspring.

# Arguments
- `ll`: the full line-list `DataFrame` with columns including
  `onset_date`, `patient_id`, and `source_case`.
- `obs_date`: the real-time observation cut-off `Date`. Cases with onset
  after this date are dropped.
"""
function filter_realtime(ll, obs_date::Date)
    _filter_linelist(ll, ll.onset_date .<= obs_date)
end

"""
$(TYPEDSIGNATURES)

**Counterfactual retrospective** view: trim to cases known to have been
infected by `obs_date` (sourced cases by their exposure upper bound,
index cases by onset).
Source attributions outside the retained set are dropped and `Z` is
rebuilt.
Relies on eventual exposure attribution and is *not* realisable in real
time; use as a comparator for the corrected real-time fit only.

# Arguments
- `ll`: the full line-list `DataFrame` with columns including
  `onset_date`, `exposure_upper`, `patient_id`, and `source_case`.
- `obs_date`: the cut-off `Date`. Sourced cases are retained when
  `exposure_upper ≤ obs_date`; index cases when `onset_date ≤ obs_date`.
"""
function filter_by_exposure(ll, obs_date::Date)
    mask = map(eachrow(ll)) do r
        ismissing(r.exposure_upper) ? r.onset_date <= obs_date :
        r.exposure_upper <= obs_date
    end
    return _filter_linelist(ll, mask)
end

function _filter_linelist(ll, mask::AbstractVector{Bool})
    sub = ll[mask, :]
    kept_ids = Set(sub.patient_id)
    src_raw = passmissing(_parse_source).(sub.source_case)
    sub.source_case = [ismissing(s) || !(string(s) in kept_ids) ?
                       missing : string(s) for s in src_raw]
    id_to_row = Dict(r.patient_id => i for (i, r) in enumerate(eachrow(sub)))
    Z = zeros(Int, nrow(sub))
    for r in eachrow(sub)
        ismissing(r.source_case) && continue
        haskey(id_to_row, r.source_case) || continue
        Z[id_to_row[r.source_case]] += 1
    end
    sub.Z = Z
    return sub
end

#
# Convolved-delay distribution `δ + Inc(sec)` and its cdf
#
# cdf(d, Δ) = P(δ + Inc(sec) ≤ Δ)
#           = ∫ cdf(d.inc, Δ − δ) · pdf(d.δ, δ) dδ
#
# Implemented through Integrals.jl's `IntegralProblem` with a fixed
# reference domain `(-1, 1)` and the parametric integration bounds
# (`mean(d.δ) ± K·std(d.δ)`) passed via the `params` field. Change of
# variable `δ = μ + halfwidth · u` lives inside the integrand. This
# keeps the IntegralsMooncakeExt rrule on the supported `p`-tangent
# path (bounds-tangents in Integrals.jl are Zygote-only), while still
# routing through Integrals.jl so the choice of GL nodes/weights and
# any future algorithm swap stays a library concern.

const _CONVOLVED_DELAYS_ALG = GaussLegendre(; n = 80)
const _CONVOLVED_DELAYS_K = 20

function _convolved_delays_integrand(u::Real, p)
    δ = p.μ + p.halfwidth * u
    return [(x_i - δ) > 0 ?
            cdf(p.inc_dist, x_i - δ) * pdf(p.δ_dist, δ) :
            zero(p.halfwidth)
            for x_i in p.x]
end

"""
$(TYPEDSIGNATURES)

Distribution of the convolved chain delay `δ + Inc(sec)`, with
`Inc(sec) ~ inc` and `δ ~ δ`. Used as the joint right-truncation
distribution for sourced offspring at a real-time cut-off.

Two `cdf` methods are provided. The scalar form `cdf(d, x::Real)` is
the standard `Distributions` interface. The vector form
`cdf(d, xs::AbstractVector)` evaluates the CDF at every `x ∈ xs` in a
single GaussLegendre solve. Prefer the vector form when evaluating at
many points: `cdf.(d, xs)` would otherwise trigger one solve per
element.

# Fields
- `inc`: incubation period distribution for the secondary case.
- `δ`: per-pair transmission timing distribution.
"""
struct ConvolvedDelays{I, D} <: ContinuousUnivariateDistribution
    inc::I
    δ::D
end

Distributions.cdf(d::ConvolvedDelays, x::Real) = cdf(d, [x])[1]
function Distributions.cdf(d::ConvolvedDelays, x::AbstractVector)
    μ = mean(d.δ)
    halfwidth = _CONVOLVED_DELAYS_K * std(d.δ)
    params = (; μ, halfwidth, x, inc_dist = d.inc, δ_dist = d.δ)
    prob = IntegralProblem(_convolved_delays_integrand, (-1.0, 1.0), params)
    return halfwidth .* solve(prob, _CONVOLVED_DELAYS_ALG).u
end

# Closed-form specialisations: defer to `Distributions.convolve` for
# same-family pairs that have an analytic sum, skipping the integral.
function Distributions.cdf(d::ConvolvedDelays{<:Normal, <:Normal}, x::Real)
    return cdf(Distributions.convolve(d.inc, d.δ), x)
end
function Distributions.cdf(d::ConvolvedDelays{<:Normal, <:Normal},
        x::AbstractVector)
    return cdf.(Distributions.convolve(d.inc, d.δ), x)
end

Distributions.minimum(::ConvolvedDelays) = -Inf
Distributions.maximum(::ConvolvedDelays) = Inf

"""
$(TYPEDSIGNATURES)

Evaluate the time-varying reproduction number at each per-case onset
offset `T_d[i]` using piecewise-constant interpolation against `edges`
and the per-draw log-R vector `logR_d`.
Returns a `Vector{Float64}` of length `length(T_d)`.
"""
function _rates_at_onsets(T_d, edges, logR_d)
    return [exp(log_R_at(T_d[i], edges, logR_d))
            for i in eachindex(T_d)]
end

"""
$(TYPEDSIGNATURES)

Per-source probability that the offspring chain `δ + Inc(sec)` has
completed by elapsed time `Δ` (vector of per-source horizons measured
in days from the source's onset). Used by
[`truncation_model`](@ref)'s real-time denominator.
"""
_chain_completion(inc, δd, Δ) = cdf(ConvolvedDelays(inc, δd), Δ)

function _pipeline_integrand(u::Real, p)
    δ = p.mid + p.halfwidth * u
    return (p.Δ_p - δ) > 0 ?
           ccdf(p.inc_dist, p.Δ_p - δ) * pdf(p.δ_dist, δ) :
           pdf(p.δ_dist, δ)
end

"""
$(TYPEDSIGNATURES)

Probability that an offspring of a source case is currently in
incubation at `Δ_p` days after the source's onset, given the
counterfactual that transmission from that source stops `Δ_q` days
after its onset:

```
P(δ ≤ Δ_q ∧ δ + Inc > Δ_p)
    = ∫_{-∞}^{Δ_q} pdf(δ_dist, δ) · ccdf(inc, Δ_p − δ) dδ
```

(where `ccdf(inc, x) = 1` for `x ≤ 0`). Evaluated by GaussLegendre
quadrature over the truncated support
`(mean(δ_dist) − K·std(δ_dist), Δ_q)`. The natural-chain limit
(no intervention) falls out as `Δ_q = +Inf`, where the upper bound is
clipped to `mean(δ_dist) + K·std(δ_dist)` and the result reduces to
`1 − cdf(ConvolvedDelays(inc, δ_dist), Δ_p)`. The exclusion sentinel
`Δ_q = -Inf` returns exactly `0`.

# Arguments
- `inc`: incubation period distribution for the secondary case.
- `δd`: per-pair transmission-timing distribution `δ`.
- `Δ_q`: per-source intervention offset (days from source onset).
- `Δ_p`: per-source observation offset (days from source onset).
"""
function _pipeline_probability(inc, δd, Δ_q::Real, Δ_p::Real)
    μ = mean(δd)
    σ = std(δd)
    a = μ - _CONVOLVED_DELAYS_K * σ
    b = min(Δ_q, μ + _CONVOLVED_DELAYS_K * σ)
    b <= a && return zero(Δ_p)
    halfwidth = (b - a) / 2
    mid = (a + b) / 2
    params = (; mid, halfwidth, Δ_p, inc_dist = inc, δ_dist = δd)
    prob = IntegralProblem(_pipeline_integrand, (-1.0, 1.0), params)
    return halfwidth * solve(prob, _CONVOLVED_DELAYS_ALG).u
end

# Resolve `intervention_time` into per-source offsets `Δ_q[i]` measured
# from each source's own onset time `T_d[i]`, in units of days from `t0`.
# The `:natural` sentinel resolves to all-`+Inf` offsets (the no-
# intervention natural-chain counterfactual).
function _intervention_offsets(intervention_time, obs_offset, T_d, t0, N)
    if intervention_time === :natural
        return fill(Inf, N)
    elseif intervention_time === nothing
        return [obs_offset - T_d[i] for i in 1:N]
    elseif intervention_time isa Date
        Δ_int = Float64(Dates.value(intervention_time - t0))
        return [Δ_int - T_d[i] for i in 1:N]
    elseif intervention_time isa AbstractVector
        length(intervention_time) == N || throw(ArgumentError(
            "intervention_time vector length $(length(intervention_time)) " *
            "does not match number of observed sources $(N)"))
        return map(1:N) do i
            it = intervention_time[i]
            it === nothing && return -Inf
            it isa Date || throw(ArgumentError(
                "intervention_time vector must contain Date or nothing"))
            Float64(Dates.value(it - t0)) - T_d[i]
        end
    else
        throw(ArgumentError(
            "intervention_time must be nothing, a Date, " *
            "or a Vector of Date/Nothing"))
    end
end

function _predict_future_onsets(model, chn, post, d;
        obs_time::Date, t0::Date,
        intervention_time = nothing,
        rng = Random.MersenneTwister(2026))
    edges = prepare_rt_edges(d.t0; obs_time = obs_time)
    obs_offset = eltype(edges)(Dates.value(obs_time - d.t0))

    (; μ_inc, σ_inc, μ_δ, σ_δ, k) = post
    log_R = post.log_R_chain
    T_on = vector_chain(chn, :T_onset)

    inc_sub = model.defaults.incubation
    δ_sub = model.defaults.transmission
    cases_sub = model.defaults.cases

    n_draws = length(k)
    N = d.N
    Z_obs = d.Zobs
    future_total = Vector{Int}(undef, n_draws)

    for d_idx in 1:n_draws
        inc = DynamicPPL.fix(inc_sub,
            (; μ_inc = μ_inc[d_idx], σ_inc = σ_inc[d_idx]))().dist
        δd = DynamicPPL.fix(δ_sub,
            (; μ_δ = μ_δ[d_idx], σ_δ = σ_δ[d_idx]))().dist
        k_d = k[d_idx]
        logR_d = [log_R[b][d_idx] for b in eachindex(log_R)]
        T_d = [T_on[i][d_idx] for i in 1:N]

        Δ_p = [obs_offset - T_d[i] for i in 1:N]
        Δ_q = _intervention_offsets(intervention_time, obs_offset,
            T_d, d.t0, N)
        p_vec = _chain_completion(inc, δd, Δ_p)
        R_vec = _rates_at_onsets(T_d, edges, logR_d)

        zsum = 0
        for i in 1:N
            R_i = R_vec[i]
            p_i = p_vec[i]
            prob = _pipeline_probability(inc, δd, Δ_q[i], Δ_p[i])
            zsum += posterior_predictive(cases_sub, rng,
                k_d, Z_obs[i], R_i, p_i, prob)
        end
        future_total[d_idx] = zsum
    end

    return (; future_samples = future_total, n_obs = N)
end

"""
$(TYPEDSIGNATURES)

Count line-list cases with onset strictly after `obs_time`. The
realised analogue of the `future_samples` field returned by
[`predict_controlled_outbreak`](@ref) and
[`predict_natural_chain_outbreak`](@ref); kept as a separate function
so the predictors stay pure functions of the fit.

# Arguments
- `ll`: full line-list `DataFrame` with an `onset_date` column.
- `obs_time`: real-time cut-off `Date`. Cases whose onset is strictly
  after this date are counted.
"""
function realised_future_count(ll, obs_time::Date)
    sum(skipmissing(ll.onset_date) .> obs_time)
end

"""
$(TYPEDSIGNATURES)

Strict controlled-outbreak counterfactual: predict future onsets
assuming all transmission stops at `intervention_time` (default
`obs_time`).
Only people already infected by `intervention_time` (transmission event
before the cut-off, chain not yet symptomatic) contribute.
For each posterior draw and observed source `i`, the contribution is
`Poisson(λ_i · pipeline_i)` where
`pipeline_i = P(δ ≤ Δ_q[i] ∧ δ + Inc > Δ_p[i])`
is the joint probability that transmission from source `i` happened
by `Δ_q[i] = intervention_time_i − T_onset[i]` *and* the offspring is
still in incubation at the observation horizon
`Δ_p[i] = obs_time − T_onset[i]`,
`p_i = cdf(ConvolvedDelays(inc, δ), Δ_p[i])` is the probability the
full chain has completed by `obs_time`, and
`λ_i | Z_obs[i], k, R_i, p_i ~ Gamma(k + Z_obs[i], R_i / (k + R_i · p_i))`.
`pipeline_i` is evaluated via [`_pipeline_probability`](@ref), which
truncates `δ` at `Δ_q[i]` before convolving with `Inc` so the formula
is correct whether the intervention precedes, coincides with, or
follows the observation cut-off.
Incubation and transmission timing distributions per draw are
recovered from the fitted `model`'s submodels via `DynamicPPL.fix`.

The prediction is a pure function of the fit (the fitted `model`, its
chain, posterior summary, and the same `d` the model was fit on).
[`realised_future_count`](@ref) returns the corresponding realised
count from the full line list as a separate evaluation step.

# Arguments
- `model`: the `DynamicPPL.Model` that produced `chn`, carrying the
  incubation, transmission, and observation submodels under
  `model.defaults`.
- `chn`: chain from a real-time joint fit.
- `post`: posterior summary returned by `summarise(chn)`.
- `d`: structured line-list data the model was fit on (returned by
  `build_data`).

# Keyword Arguments
- `obs_time`: cut-off `Date`.
- `t0`: time origin `Date`.
- `intervention_time`: when transmission stops. Three forms are
  accepted:
  - `nothing` (default): intervention coincides with `obs_time` for
    every source, equivalent to passing `obs_time` as a scalar `Date`.
  - A scalar `Date`: applies uniformly to every source.
    `Δ_q[i] = intervention_time − T_onset[i]`. Sources whose onset is
    after `intervention_time` get a negative `Δ_q[i]`, so
    `q_i = cdf(δ, Δ_q[i])` evaluates in the left tail of the
    transmission-timing distribution.
  - A `Vector{Union{Date, Nothing}}` of length `N` (the number of
    observed sources in `d`): per-source intervention dates. A
    `nothing` element means the source contributes zero future onsets
    (resolved as `q_i = cdf(δ, -Inf) = 0`), useful when the caller
    wants to encode "only these sources had ongoing transmission" as
    policy in the calling code rather than in the package.
- `rng`: RNG used for posterior-predictive draws.

# Example

Per-source intervention dates for case-by-case isolation:

```julia
isolation = [c.onset_date + Day(1) for c in eachrow(ll_rt)]
predict_controlled_outbreak(model, chn, post, d_rt;
    obs_time = obs_date, t0 = t0, intervention_time = isolation)
```
"""
function predict_controlled_outbreak(model, chn, post, d;
        obs_time::Date, t0::Date,
        intervention_time = nothing,
        rng = Random.MersenneTwister(2026))
    return _predict_future_onsets(model, chn, post, d;
        obs_time = obs_time, t0 = t0,
        intervention_time = intervention_time, rng = rng)
end

"""
$(TYPEDSIGNATURES)

Natural-chain counterfactual: predict future onsets assuming
currently-observed sources keep transmitting at their existing rate
but no second-generation chains form from those new offspring.
Each observed source contributes `Poisson(λ_i · (1 − p_i))` future
onsets, covering both offspring already infected by `obs_time` (chain
in incubation) and offspring still to be infected by the source.
Same Gamma posterior on `λ_i` as [`predict_controlled_outbreak`](@ref);
this is the `Δ_q = +Inf` limit of that predictor, where
[`_pipeline_probability`](@ref) reduces to `1 − p_i`.
Incubation and transmission timing distributions per draw are
recovered from the fitted `model`'s submodels via `DynamicPPL.fix`.

# Arguments
- `model`: the `DynamicPPL.Model` that produced `chn`, carrying the
  incubation, transmission, and observation submodels under
  `model.defaults`.
- `chn`: chain from a real-time joint fit.
- `post`: posterior summary returned by `summarise(chn)`.
- `d`: structured line-list data the model was fit on (returned by
  `build_data`).

# Keyword Arguments
- `obs_time`: cut-off `Date`.
- `t0`: time origin `Date`.
- `rng`: RNG used for posterior-predictive draws.
"""
function predict_natural_chain_outbreak(model, chn, post, d;
        obs_time::Date, t0::Date,
        rng = Random.MersenneTwister(2026))
    return _predict_future_onsets(model, chn, post, d;
        obs_time = obs_time, t0 = t0,
        intervention_time = :natural, rng = rng)
end
