## Real-time extension utilities.
##
## Helpers to simulate a real-time analysis from a closed-out line list and
## the cluster-completeness integral that goes into the right-truncated
## likelihood.

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
# δ ~ Normal(μ_δ, σ_δ). Evaluated for every source case on every NUTS
# gradient step.
#
# Implemented as a 2-D quadrature in standardised normal coordinates
# (z₁ = log-Inc(src), z₂ = δ): the change of variables removes the
# parameter-dependent integration domain and bounds the integrand. The
# inner CDF is the LogNormal CDF of the secondary's incubation.
#
# AD path: Mooncake reverse mode via DifferentiationInterface.jl, which
# consumes the ChainRulesCore `rrule` that `Integrals.jl` defines for
# `__solvebp` (Enzyme reverse mode has no equivalent EnzymeRule and so
# is incompatible with the Integrals.jl solve machinery).

# Outer integration box in standardised coordinates. ±8σ captures the
# integrand to better than 1e-10 in the model regime.
const _F_CLUSTER_BOX = ([-8.0, -8.0], [8.0, 8.0])

# Default algorithm: adaptive 2-D HCubature, tight tolerance. Override
# via the `fcluster_alg` kwarg to `joint_model` and `analyse` for
# diagnostics or to trade accuracy for speed.
const DEFAULT_FCLUSTER_ALG = HCubatureJL()

function _f_cluster_integrand(z::AbstractVector, p)
    z1, z2 = z[1], z[2]
    inc_src = exp(p.μ_inc + p.σ_inc * z1)
    δ       = p.μ_δ + p.σ_δ * z2
    s       = p.t - inc_src - δ
    ϕ1 = exp(-z1 * z1 / 2) / sqrt(2π)
    ϕ2 = exp(-z2 * z2 / 2) / sqrt(2π)
    inc_cdf = s > 0 ? cdf(LogNormal(p.μ_inc, p.σ_inc), s) : zero(s)
    return inc_cdf * ϕ1 * ϕ2
end

"""
    F_cluster(t, μ_inc, σ_inc, μ_δ, σ_δ;
              alg = HCubatureJL(), reltol = 1e-8, abstol = 1e-10)

Probability that the full source-to-secondary chain
`Inc(src) + δ + Inc(sec)` is no greater than `t`, with
`Inc(src), Inc(sec) ~ LogNormal(μ_inc, σ_inc)` (i.i.d.) and
`δ ~ Normal(μ_δ, σ_δ)`.

`alg` is any `Integrals.AbstractIntegralAlgorithm`; the default is
`HCubatureJL()`. Use a different algorithm (e.g. `QuadratureRule` with
fixed Gauss-Hermite nodes) to trade adaptive accuracy for evaluation
speed during diagnostics.
"""
function F_cluster(t::Real, μ_inc::Real, σ_inc::Real, μ_δ::Real, σ_δ::Real;
                   alg = DEFAULT_FCLUSTER_ALG,
                   reltol::Real = 1e-8, abstol::Real = 1e-10)
    t <= 0 && return zero(float(t))
    p = (; t, μ_inc, σ_inc, μ_δ, σ_δ)
    prob = IntegralProblem(_f_cluster_integrand, _F_CLUSTER_BOX, p)
    sol = solve(prob, alg; reltol = reltol, abstol = abstol)
    return sol.u
end

"""
    F_cluster_vec(θ; alg = HCubatureJL())

Vector-input wrapper of [`F_cluster`](@ref) for AD testing.
`θ = [t, μ_inc, σ_inc, μ_δ, σ_δ]`.
"""
F_cluster_vec(θ; alg = DEFAULT_FCLUSTER_ALG) =
    F_cluster(θ[1], θ[2], θ[3], θ[4], θ[5]; alg = alg)
