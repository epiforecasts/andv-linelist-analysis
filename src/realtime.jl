## Real-time extension utilities.
##
## Helpers to simulate a real-time analysis from a closed-out line list and
## the offspring-completeness integral that goes into the right-truncated
## likelihood.

using FastGaussQuadrature: gausshermite
using Integrals: QuadratureRule
using SpecialFunctions: erfc

"""
    filter_realtime(ll, obs_date)

Return a copy of the line-list trimmed to cases whose onset is on or
before `obs_date`, i.e. the cases an analyst would know about at that
cut-off. Source attributions that point outside the retained set are
dropped (the attributed case becomes an apparent index) and the `Z`
column is rebuilt to count only retained offspring.

The returned line list is suitable for feeding back into
[`build_data`](@ref) with `obs_time = obs_date`.
"""
function filter_realtime(ll, obs_date::Date)
    return _filter_linelist(ll, ll.onset_date .<= obs_date)
end

"""
    filter_by_exposure(ll, obs_date)

Return a copy of the line-list trimmed to cases known to have been
infected by `obs_date`. This is a **counterfactual retrospective** view:
for sourced cases we keep those whose exposure upper bound is on or
before `obs_date`; for index cases (no recorded exposure window) we
keep those whose onset is on or before `obs_date`. Source attributions
that point outside the retained set are dropped and `Z` is rebuilt.

This filter relies on information not available in real time —
specifically, the eventual exposure attribution for cases still in
incubation at `obs_date`. It is intended as a comparator for the
corrected real-time fit, not as a real analysis path.
"""
function filter_by_exposure(ll, obs_date::Date)
    mask = map(eachrow(ll)) do r
        if ismissing(r.exposure_upper)
            r.onset_date <= obs_date
        else
            r.exposure_upper <= obs_date
        end
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

# ---------------------------------------------------------------------------
# Offspring-completeness integral
# ---------------------------------------------------------------------------
#
# F_offspring(Δ; μ_inc, σ_inc, μ_δ, σ_δ) = P(δ + Inc(sec) ≤ Δ)
#
# where Inc(sec) ~ LogNormal(μ_inc, σ_inc) and δ ~ Normal(μ_δ, σ_δ).
# Evaluated for every observed source case on every NUTS gradient step,
# with Δ = obs_time − T_onset[src]. The source's own incubation is *not*
# marginalised over here: in `joint_model` `Inc(src) = T_onset[src] −
# T_inf[src]` is a sampled latent already scored against `inc_dist`, so
# conditioning on it leaves only `δ` and `Inc(sec)` free for the
# offspring-completeness probability.
#
# Implemented as a 1-D quadrature in the standardised normal coordinate
# z attached to the Normal δ axis,
#
#   δ = μ_δ + σ_δ · z,   z ~ N(0, 1),
#
# with the LogNormal Inc(sec) CDF in closed form:
#
#   F_offspring(Δ) = ∫ F_inc(Δ − μ_δ − σ_δ z) ϕ(z) dz,
#
# where `F_inc` is the LogNormal(μ_inc, σ_inc) CDF (returning 0 when its
# argument is non-positive). Integrating over the δ axis (not the
# Inc(sec) axis) keeps the integrand smooth: `F_inc(t − δ)` is C∞ in `δ`
# wherever the argument is positive, so Gauss-Hermite converges
# essentially exactly at modest node counts. The mirror formulation —
# integrating over the Inc(sec) log-axis with the Normal δ CDF in
# closed form — looks superficially symmetric but has a step-like
# integrand (Φ composed with exp) and converges very slowly.
#
# A single vector-valued solve evaluates the integral at every requested
# Δ simultaneously, so one quadrature replaces N scalar solves in the
# model loop.
#
# Solver: `Integrals.QuadratureRule` with pre-computed Gauss-Hermite
# nodes. A fixed-rule solver is required for Mooncake reverse-mode AD:
# `HCubatureJL`'s adaptive heap mutates internal arrays in ways that
# trip an unhandled LLVM intrinsic (`sub_ptr`) inside Mooncake, even
# when the integrand itself is differentiation-friendly. Gauss-Hermite
# is the natural rule because the standardised coordinate already
# carries a Gaussian weight.

# Pre-compute Gauss-Hermite nodes once at module load. Returning const
# arrays from the quadrature `q(n)` callback avoids pulling
# `FastGaussQuadrature.gausshermite` (and its BLAS eigensolver) into the
# differentiated path.
const _GH_N = 14
const _GH_NODES, _GH_WEIGHTS = let
    nodes_1d, weights_1d = gausshermite(_GH_N)
    # Change of variables to a standard normal: δ = μ_δ + σ_δ z,
    # so absorb sqrt(2) into the nodes and divide weights by sqrt(π).
    n = sqrt(2) .* nodes_1d
    w = weights_1d ./ sqrt(π)
    # Wrap each node in a 1-element vector so the QuadratureRule callback
    # matches the `AbstractVector` signature Integrals.jl expects.
    nodes   = [Float64[x] for x in n]
    weights = collect(w)
    nodes, weights
end
_gh_q(::Int) = (_GH_NODES, _GH_WEIGHTS)

# `evalrule` in Integrals.jl rescales nodes by (ub-lb)/2 and shifts by
# (ub+lb)/2; using [-1, 1] means scale = 1, shift = 0, so the rule
# evaluates the integrand at the pre-computed nodes unchanged.
const _F_OFFSPRING_DOMAIN = ([-1.0], [1.0])
const _F_OFFSPRING_ALG    = QuadratureRule(_gh_q; n = _GH_N)

# Closed-form LogNormal CDF, avoids constructing a `LogNormal` per node.
# cdf(LogNormal(μ, σ), x) = 0.5 * erfc(-(log(x) - μ) / (σ * √2))
@inline _lognormal_cdf(x, μ, σ) =
    x > 0 ? oftype(x, 0.5) * erfc(-(log(x) - μ) / (σ * sqrt(oftype(x, 2)))) : zero(x)

# Out-of-place vector-output integrand. At one z point compute the
# realised δ, then evaluate the LogNormal CDF of Inc(sec) at every
# `Δ − δ` requested in `p.ts`; returns a fresh length-N vector.
function _f_offspring_integrand(z::AbstractVector, p)
    z1 = z[1]
    δ  = p.μ_δ + p.σ_δ * z1
    # ϕ(z) is already absorbed into the Gauss-Hermite weights.
    return [_lognormal_cdf(p.ts[i] - δ, p.μ_inc, p.σ_inc)
            for i in eachindex(p.ts)]
end

"""
    F_offspring(ts, inc_dist::LogNormal, δ_dist::Normal; alg = _F_OFFSPRING_ALG)

Probability that the offspring's transmission-plus-incubation chain
`δ + Inc(sec)` is no greater than `Δ`, evaluated at every element of
`ts` against the same population distributions via a single 1-D
quadrature.

`Inc(sec) ~ inc_dist` and `δ ~ δ_dist`. Returns a vector of length
`length(ts)`.

The integral is the conditional offspring-completeness for an observed
source whose onset time `T_onset[src]` is already pinned by the model:
the caller supplies `Δ = obs_time − T_onset[src]`, and the source's own
incubation is *not* marginalised over here (it is a sampled latent
scored elsewhere in `joint_model`).

`alg` is any `Integrals.AbstractIntegralAlgorithm`. Default is a 14-point
Gauss-Hermite `QuadratureRule` with nodes precomputed at module load.
Pass e.g. `HCubatureJL()` to swap in an adaptive solver for accuracy
diagnostics (note: HCubatureJL is not compatible with Mooncake
reverse-mode AD).
"""
function F_offspring(ts::AbstractVector,
                     inc_dist::LogNormal, δ_dist::Normal;
                     alg = _F_OFFSPRING_ALG)
    μ_inc, σ_inc = inc_dist.μ, inc_dist.σ
    μ_δ,   σ_δ   = δ_dist.μ,   δ_dist.σ
    p    = (; ts, μ_inc, σ_inc, μ_δ, σ_δ)
    prob = IntegralProblem(_f_offspring_integrand, _F_OFFSPRING_DOMAIN, p)
    sol  = solve(prob, alg)
    return sol.u
end

"""
    F_offspring(t::Real, inc_dist, δ_dist; kw...)

Scalar convenience: a single-`t` query that returns the value rather
than a length-1 vector.
"""
F_offspring(t::Real, inc_dist::LogNormal, δ_dist::Normal; kw...) =
    F_offspring([float(t)], inc_dist, δ_dist; kw...)[1]

"""
    F_offspring_vec(θ; alg = _F_OFFSPRING_ALG)

Vector-input wrapper of scalar [`F_offspring`](@ref) for AD testing.
`θ = [t, μ_inc, σ_inc, μ_δ, σ_δ]`; the distributions are reconstructed
internally so callers (e.g. DifferentiationInterface.jl) keep working
with a flat real-valued input.
"""
F_offspring_vec(θ; alg = _F_OFFSPRING_ALG) =
    F_offspring(θ[1], LogNormal(θ[2], θ[3]), Normal(θ[4], θ[5]); alg = alg)
