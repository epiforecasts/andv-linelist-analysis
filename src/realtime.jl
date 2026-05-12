## Real-time extension utilities.
##
## Helpers to simulate a real-time analysis from a closed-out line list and
## the cluster-completeness integral that goes into the right-truncated
## likelihood.

"""
    filter_realtime(ll, obs_date)

Return a copy of the line-list trimmed to cases whose onset is on or before
`obs_date`, i.e. the cases an analyst would know about at that cut-off.
Source attributions that point outside the retained set are dropped (the
attributed case becomes an apparent index), and the `Z` column is rebuilt
to count only retained offspring.

The returned line list is suitable for feeding back into [`build_data`](@ref)
with `obs_time = obs_date`.
"""
function filter_realtime(ll, obs_date::Date)
    keep = ll.onset_date .<= obs_date
    sub  = ll[keep, :]
    kept_ids = Set(sub.patient_id)
    src_raw  = passmissing(_parse_source).(sub.source_case)
    sub.source_case = [ismissing(s) || !(string(s) in kept_ids) ?
                       missing : string(s) for s in src_raw]
    # Rebuild Z to count retained offspring only.
    id_to_row = Dict(r.patient_id => i for (i, r) in enumerate(eachrow(sub)))
    Z = zeros(Int, nrow(sub))
    for r in eachrow(sub)
        ismissing(r.source_case) && continue
        src = r.source_case
        haskey(id_to_row, src) || continue
        Z[id_to_row[src]] += 1
    end
    sub.Z = Z
    return sub
end
##
## A real-time fit conditions every case on a per-case observation cut-off
## `obs_time[i]`. Two corrections enter the likelihood:
##
##   1. Right-truncation of Inc and δ — handled in `model.jl` by subtracting
##      `logcdf(...)` from each per-case contribution.
##   2. Cluster-completeness on the NB offspring count — the source's true
##      offspring chain `Inc(src) + δ + Inc(sec)` must fit within
##      `Δ = obs_time[src] - T_inf[src]` for the secondary to be observed.
##
## `F_cluster(Δ; μ_inc, σ_inc, μ_δ, σ_δ)` is the population-level probability
## that the chain fits in Δ. Two implementations live here:
##
##   - `F_cluster`         — tensor Gauss-Hermite over standardised (z₁, z₂).
##                           Static loops, no allocations, no data-dependent
##                           branches in the differentiated path → clean under
##                           Enzyme via DifferentiationInterface.jl. This is
##                           what the joint model calls on every NUTS step.
##   - `F_cluster_quadrature` — Integrals.jl + HCubatureJL, kept as a
##                              high-accuracy reference for unit-testing the
##                              fixed-node version.

# Gauss-Hermite nodes/weights for the physicist Hermite weight exp(-z²).
# To integrate g(z) ϕ(z) over z, evaluate at √2·node with weight / √π
# (change of variables from physicist to probabilist Hermite).
const _GH_N = 40
const _GH_NODES, _GH_WEIGHTS = let
    nodes, weights = gausshermite(_GH_N)
    sqrt(2) .* nodes, weights ./ sqrt(π)
end

# LogNormal CDF written explicitly via erf — keeps the integrand allocation-
# free and removes a struct construction from the Enzyme tape.
_lognormal_cdf(s, μ, σ) =
    s > 0 ? (one(s) + erf((log(s) - μ) / (σ * sqrt(2)))) / 2 : zero(s)

"""
    F_cluster(t, μ_inc, σ_inc, μ_δ, σ_δ)

Probability that the full source-to-secondary chain
`Inc(src) + δ + Inc(sec)` is no greater than `t`, with
`Inc(src), Inc(sec) ~ LogNormal(μ_inc, σ_inc)` (i.i.d.) and
`δ ~ Normal(μ_δ, σ_δ)`.

Tensor 20-node Gauss-Hermite quadrature in standardised normal coordinates.
Smooth in all parameters; differentiable under Enzyme via
`DifferentiationInterface.jl`.
"""
function F_cluster(t::Real, μ_inc::Real, σ_inc::Real, μ_δ::Real, σ_δ::Real)
    # Promote to a common element type so the AD tape sees a single type
    # throughout the loop.
    T = float(promote_type(typeof(t), typeof(μ_inc), typeof(σ_inc),
                           typeof(μ_δ), typeof(σ_δ)))
    t <= 0 && return zero(T)
    acc = zero(T)
    @inbounds for i in eachindex(_GH_NODES)
        z1 = _GH_NODES[i]
        w1 = _GH_WEIGHTS[i]
        inc_src = exp(μ_inc + σ_inc * z1)
        for j in eachindex(_GH_NODES)
            z2 = _GH_NODES[j]
            w2 = _GH_WEIGHTS[j]
            δ  = μ_δ + σ_δ * z2
            s  = t - inc_src - δ
            acc += w1 * w2 * _lognormal_cdf(s, μ_inc, σ_inc)
        end
    end
    return acc
end

"""
    F_cluster_quadrature(t, μ_inc, σ_inc, μ_δ, σ_δ; reltol = 1e-9, abstol = 1e-12)

Reference implementation of [`F_cluster`](@ref) using `Integrals.jl` with
`HCubatureJL` as the back-end. Two orders of magnitude slower than the
fixed-node default but adaptive to higher accuracy; kept for unit tests
that pin the fixed-node result against an independent solver.
"""
function F_cluster_quadrature(t::Real, μ_inc::Real, σ_inc::Real,
                              μ_δ::Real, σ_δ::Real;
                              reltol::Real = 1e-9, abstol::Real = 1e-12,
                              z_bound::Real = 8.0)
    t <= 0 && return zero(float(t))
    p = (; t, μ_inc, σ_inc, μ_δ, σ_δ)
    prob = IntegralProblem(_f_cluster_integrand_z,
                           ([-z_bound, -z_bound], [z_bound, z_bound]), p)
    sol = solve(prob, HCubatureJL(); reltol = reltol, abstol = abstol)
    return sol.u
end

function _f_cluster_integrand_z(z::AbstractVector, p)
    z1, z2 = z[1], z[2]
    inc_src = exp(p.μ_inc + p.σ_inc * z1)
    δ       = p.μ_δ + p.σ_δ * z2
    s       = p.t - inc_src - δ
    ϕ1 = exp(-z1 * z1 / 2) / sqrt(2π)
    ϕ2 = exp(-z2 * z2 / 2) / sqrt(2π)
    return _lognormal_cdf(s, p.μ_inc, p.σ_inc) * ϕ1 * ϕ2
end

# AD-test wrappers: DifferentiationInterface expects f(x::AbstractVector).
"""
    F_cluster_vec(θ) -> Real

Vector-input wrapper of [`F_cluster`](@ref) for AD testing.
`θ = [t, μ_inc, σ_inc, μ_δ, σ_δ]`.
"""
F_cluster_vec(θ) = F_cluster(θ[1], θ[2], θ[3], θ[4], θ[5])

"""
    F_cluster_quadrature_vec(θ) -> Real

Vector-input wrapper of [`F_cluster_quadrature`](@ref) for AD testing.
"""
F_cluster_quadrature_vec(θ) = F_cluster_quadrature(θ[1], θ[2], θ[3], θ[4], θ[5])
