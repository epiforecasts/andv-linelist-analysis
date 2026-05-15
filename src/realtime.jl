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
Distributions.minimum(::ConvolvedDelays) = -Inf
Distributions.maximum(::ConvolvedDelays) = Inf

"""
$(TYPEDSIGNATURES)

Controlled-outbreak counterfactual: for each posterior draw, sum the
predicted future onsets across observed sources, assuming no further
transmission after `obs_time`. Each observed source contributes a
`Poisson(λ_i · (1 − p_i))` term where `λ_i | Z_obs[i], k, R_i, p_i`
follows the conjugate Gamma posterior of the NB-binomial-thinning
model: `Gamma(k + Z_obs[i], R_i / (k + R_i · p_i))` (scale form),
with `p_i = cdf(ConvolvedDelays(inc, δ), obs_time − T_onset[i])`. Conditioning on each
source's observed offspring count `Z_obs[i]` sharpens the prediction
relative to the naive marginal NB(k, R_i (1 − p_i)).

Also returns the count of cases with onset strictly after `obs_time`
in the supplied (full) line list, as a comparator for the
counterfactual: the realised count lies below the predicted band if
control were achieved at the cut-off, above if transmission continued.

# Arguments
- `chn`: chain from a real-time `analyse(; obs_time, ...)` fit.
- `post`: posterior summary returned by `summarise(chn)` (or
  the second return value of `analyse`).
- `ll`: full line-list `DataFrame` (not the obs_time-filtered subset).
- `obs_time`: cut-off `Date` used to fit `chn`.
- `t0`: time origin `Date` used to fit `chn`.

# Keyword Arguments
- `inc_dist`: two-argument constructor for the incubation period
  distribution, called as `inc_dist(μ_inc, σ_inc)` per posterior draw.
  Defaults to `LogNormal`; override to match an alternative incubation
  parameterisation in the joint model.
- `delta_dist`: two-argument constructor for the transmission timing
  distribution, called as `delta_dist(μ_δ, σ_δ)`. Defaults to `Normal`.
- `rng`: RNG used for posterior-predictive draws.
"""
function predict_controlled_outbreak(chn, post, ll,
        obs_time::Date, t0::Date;
        inc_dist = LogNormal,
        delta_dist = Normal,
        rng = Random.MersenneTwister(2026))
    ll_rt = filter_realtime(ll, obs_time)
    d = build_data(ll_rt; obs_time = obs_time, t0 = t0)
    edges = bin_edges_day(d.t0)
    T = eltype(edges)
    obs_offset = T(Dates.value(obs_time - d.t0))
    edges = edges[edges .<= obs_offset]
    if isempty(edges) || edges[end] < obs_offset
        push!(edges, obs_offset)
    end

    (; μ_inc, σ_inc, μ_δ, σ_δ, k) = post
    log_R = post.log_R_chain
    T_on = vector_chain(chn, :T_onset)

    n_draws = length(k)
    N = d.N
    Z_obs = d.Zobs
    future_total = Vector{Int}(undef, n_draws)

    for d_idx in 1:n_draws
        inc = inc_dist(μ_inc[d_idx], σ_inc[d_idx])
        δd = delta_dist(μ_δ[d_idx], σ_δ[d_idx])
        k_d = k[d_idx]
        logR_d = [log_R[b][d_idx] for b in eachindex(log_R)]
        T_d = [T_on[i][d_idx] for i in 1:N]

        Δ = [obs_offset - T_d[i] for i in 1:N]
        p_vec = cdf(ConvolvedDelays(inc, δd), Δ)

        zsum = 0
        for i in 1:N
            # Clamp tightly enough that exp(log_R) stays well within
            # Int64 range under the downstream Poisson sampler. Divergent
            # draws (Mode B) can otherwise push R_i past 1e18 and overflow.
            R_i = exp(clamp(log_R_at(T_d[i], edges, logR_d),
                -T(20), T(20)))
            p_i = clamp(p_vec[i], zero(T), one(T))
            shape_post = k_d + Z_obs[i]
            scale_post = R_i / (k_d + R_i * p_i)
            λ_i = rand(rng, Gamma(shape_post, scale_post))
            zsum += rand(rng, Poisson(λ_i * (one(T) - p_i)))
        end
        future_total[d_idx] = zsum
    end

    actual_future = sum(ll.onset_date .> obs_time)
    return (; future_samples = future_total,
        actual_future,
        n_obs = N)
end
