## Real-time extension utilities.
##
## Helpers to simulate a real-time analysis from a closed-out line list and
## the cluster-completeness integral that goes into the right-truncated
## likelihood.

using FastGaussQuadrature: gausshermite
using Integrals: QuadratureRule

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
# Cluster-completeness integral
# ---------------------------------------------------------------------------
#
# F_cluster(t; μ_inc, σ_inc, μ_δ, σ_δ) = P(Inc(src) + δ + Inc(sec) ≤ t)
#
# where Inc(src), Inc(sec) ~ LogNormal(μ_inc, σ_inc) i.i.d. and
# δ ~ Normal(μ_δ, σ_δ).
# Evaluated for every source case on every NUTS gradient step.
#
# Implemented as a 2-D vector-valued quadrature in standardised normal
# coordinates (z₁ = log-Inc(src), z₂ = δ).
# The change of variables removes the parameter-dependent integration
# domain and bounds the integrand.
# The integrand returns a vector with one component per requested `t`,
# so a single solve replaces N scalar solves inside the model loop.
#
# Solver: `Integrals.QuadratureRule` with pre-computed Gauss-Hermite
# tensor nodes.
# A fixed-rule solver is required for Mooncake reverse-mode AD:
# `HCubatureJL`'s adaptive heap mutates internal arrays in ways that
# trip an unhandled LLVM intrinsic (`sub_ptr`) inside Mooncake, even
# when the integrand itself is differentiation-friendly.
# Gauss-Hermite is the natural rule here because the standardised
# coordinates already carry a Gaussian weight.

# Pre-compute Gauss-Hermite tensor nodes once at module load.
# Returning const arrays from the quadrature `q(n)` callback avoids
# pulling `FastGaussQuadrature.gausshermite` (and its BLAS eigensolver)
# into the differentiated path.
const _GH_N = 20
const _GH_NODES, _GH_WEIGHTS = let
    nodes_1d, weights_1d = gausshermite(_GH_N)
    n = sqrt(2) .* nodes_1d
    w = weights_1d ./ sqrt(π)
    nodes   = vec([Float64[x, y] for x in n, y in n])
    weights = vec([wx * wy for wx in w, wy in w])
    nodes, weights
end
_gh_tensor_q(::Int) = (_GH_NODES, _GH_WEIGHTS)

# `evalrule` in Integrals.jl rescales nodes by (ub-lb)/2 and shifts by
# (ub+lb)/2; using [-1,1]² means scale = 1, shift = 0, so the rule
# evaluates the integrand at the pre-computed nodes unchanged.
const _F_CLUSTER_DOMAIN = ([-1.0, -1.0], [1.0, 1.0])
const _F_CLUSTER_ALG    = QuadratureRule(_gh_tensor_q; n = _GH_N)

# Out-of-place vector-output integrand. At one (z₁, z₂) point evaluate
# the LogNormal CDF of the secondary's incubation at every t in p.ts,
# return a fresh length-N vector.
function _f_cluster_integrand(z::AbstractVector, p)
    z1, z2 = z[1], z[2]
    inc_src = exp(p.μ_inc + p.σ_inc * z1)
    δ       = p.μ_δ + p.σ_δ * z2
    # ϕ(z) factors are already absorbed into the Gauss-Hermite weights.
    return [let s = p.ts[i] - inc_src - δ
                s > 0 ? cdf(LogNormal(p.μ_inc, p.σ_inc), s) : zero(s)
            end
            for i in eachindex(p.ts)]
end

"""
    F_cluster(ts, μ_inc, σ_inc, μ_δ, σ_δ; alg = _F_CLUSTER_ALG)

Probability that the full source-to-secondary chain
`Inc(src) + δ + Inc(sec)` is no greater than `t`, evaluated at every
element of `ts` against the same population parameters via a single
2-D quadrature.

`Inc(src), Inc(sec) ~ LogNormal(μ_inc, σ_inc)` i.i.d., and
`δ ~ Normal(μ_δ, σ_δ)`.
Returns a vector of length `length(ts)`.

`alg` is any `Integrals.AbstractIntegralAlgorithm`.
Default is a 20×20 Gauss-Hermite tensor `QuadratureRule` with nodes
precomputed at module load.
Pass e.g. `HCubatureJL()` to swap in an adaptive solver for accuracy
diagnostics (note: HCubatureJL is not compatible with Mooncake
reverse-mode AD).
"""
function F_cluster(ts::AbstractVector,
                   μ_inc::Real, σ_inc::Real, μ_δ::Real, σ_δ::Real;
                   alg = _F_CLUSTER_ALG)
    p    = (; ts, μ_inc, σ_inc, μ_δ, σ_δ)
    prob = IntegralProblem(_f_cluster_integrand, _F_CLUSTER_DOMAIN, p)
    sol  = solve(prob, alg)
    return sol.u
end

"""
    F_cluster(t::Real, μ_inc, σ_inc, μ_δ, σ_δ; kw...)

Scalar convenience: a single-`t` query that returns the value rather
than a length-1 vector.
"""
F_cluster(t::Real, μ_inc::Real, σ_inc::Real, μ_δ::Real, σ_δ::Real; kw...) =
    F_cluster([float(t)], μ_inc, σ_inc, μ_δ, σ_δ; kw...)[1]

"""
    F_cluster_vec(θ; alg = _F_CLUSTER_ALG)

Vector-input wrapper of scalar [`F_cluster`](@ref) for AD testing.
`θ = [t, μ_inc, σ_inc, μ_δ, σ_δ]`.
"""
F_cluster_vec(θ; alg = _F_CLUSTER_ALG) =
    F_cluster(θ[1], θ[2], θ[3], θ[4], θ[5]; alg = alg)
