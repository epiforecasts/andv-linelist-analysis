## Real-time extension utilities.

import FastGaussQuadrature  # activates Integrals' FGQ extension for GaussLegendre nodes
using Integrals: GaussLegendre

"""
    filter_realtime(ll, obs_date)

Trim the line-list to cases whose onset is on or before `obs_date`
(those an analyst would know about at the cut-off). Source attributions
pointing outside the retained set are dropped and `Z` is rebuilt to count
only retained offspring.
"""
filter_realtime(ll, obs_date::Date) =
    _filter_linelist(ll, ll.onset_date .<= obs_date)

"""
    filter_by_exposure(ll, obs_date)

**Counterfactual retrospective** view: trim to cases known to have been
infected by `obs_date` (sourced cases by their exposure upper bound,
index cases by onset). Source attributions outside the retained set are
dropped and `Z` is rebuilt. Relies on eventual exposure attribution and
is *not* realisable in real time; use as a comparator for the corrected
real-time fit only.
"""
function filter_by_exposure(ll, obs_date::Date)
    mask = map(eachrow(ll)) do r
        ismissing(r.exposure_upper) ? r.onset_date <= obs_date :
                                      r.exposure_upper <= obs_date
    end
    return _filter_linelist(ll, mask)
end

function _filter_linelist(ll, mask::AbstractVector{Bool})
    sub      = ll[mask, :]
    kept_ids = Set(sub.patient_id)
    src_raw  = passmissing(_parse_source).(sub.source_case)
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
# Offspring-completeness integral
#
# F_offspring(Δ) = P(δ + Inc(sec) ≤ Δ) with δ ~ δ_dist, Inc(sec) ~ inc_dist
#                = ∫ cdf(inc_dist, Δ − δ) · pdf(δ_dist, δ) dδ
#
# Fixed-rule `GaussLegendre` rather than an adaptive solver: Mooncake
# reverse-mode AD does not survive HCubature / QuadGK's internal heap
# mutations.

function _f_offspring_integrand(δ::Real, p)
    w = pdf(p.δ_dist, δ)
    return [cdf(p.inc_dist, p.ts[i] - δ) * w for i in eachindex(p.ts)]
end

const _F_OFFSPRING_BOUNDS = (-30.0, 30.0)
const _F_OFFSPRING_ALG    = GaussLegendre(; n = 80)

"""
    F_offspring(ts, inc_dist, δ_dist;
                alg = GaussLegendre(; n = 80),
                δ_bounds = (-30.0, 30.0))

P(δ + Inc(sec) ≤ Δ) at each `Δ ∈ ts`. `Inc(sec) ~ inc_dist`,
`δ ~ δ_dist`. Returns a vector of the same length as `ts`.

The integrand calls `cdf(inc_dist, ·)` and `pdf(δ_dist, ·)` directly,
so swapping the distributional families does not require touching this
function. `δ_bounds` must contain effectively all of `δ_dist`'s mass —
the default ±30 covers any well-behaved Normal δ.
"""
function F_offspring(ts::AbstractVector, inc_dist, δ_dist;
                     alg = _F_OFFSPRING_ALG,
                     δ_bounds::Tuple{<:Real,<:Real} = _F_OFFSPRING_BOUNDS)
    p    = (; ts, inc_dist, δ_dist)
    prob = IntegralProblem(_f_offspring_integrand, δ_bounds, p)
    return solve(prob, alg).u
end

"""
    F_offspring(t::Real, inc_dist, δ_dist; kw...) -> Real

Scalar-`t` form returning the value instead of a length-1 vector.
"""
F_offspring(t::Real, inc_dist, δ_dist; kw...) =
    F_offspring([float(t)], inc_dist, δ_dist; kw...)[1]

"""
    F_offspring_vec(θ; alg = ...)

AD-test wrapper: `θ = [t, μ_inc, σ_inc, μ_δ, σ_δ]`. Reconstructs
`LogNormal` Inc and `Normal` δ from `θ` so `DifferentiationInterface`
sees a flat real-valued input.
"""
F_offspring_vec(θ; alg = _F_OFFSPRING_ALG) =
    F_offspring(θ[1], LogNormal(θ[2], θ[3]), Normal(θ[4], θ[5]); alg = alg)
