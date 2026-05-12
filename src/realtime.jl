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
    keep = ll.onset_date .<= obs_date
    sub  = ll[keep, :]
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
# δ ~ Normal(μ_δ, σ_δ). The model evaluates this on every NUTS gradient
# step, so the gradient-time path must be Enzyme reverse-mode safe via
# DifferentiationInterface.jl.
#
# The integrand is the LogNormal CDF of the secondary's incubation,
# weighted against the standardised normal density of (z₁ = log-Inc(src),
# z₂ = δ). The change of variables makes the integration domain finite
# (Gauss-Hermite nodes) and the integrand bounded.

abstract type ClusterAlg end

"""
    GaussHermite(n = 40)

Tensor-product Gauss-Hermite quadrature with `n` nodes per axis. The
loop is statically sized and uses only `+, -, *, /, exp, log, erf`, so
Enzyme reverse mode differentiates through it cleanly under
`DifferentiationInterface.jl`. This is the algorithm `joint_model` uses
on every gradient call.
"""
struct GaussHermite <: ClusterAlg
    n::Int
end
GaussHermite(; n::Int = 40) = GaussHermite(n)

"""
    HCubature(; reltol = 1e-9, abstol = 1e-12, z_bound = 8.0)

Adaptive 2-D Cubature on the (z₁, z₂) box via
`Integrals.jl` + `HCubatureJL`, integrated over `[-z_bound, z_bound]²`.
Used as a high-accuracy reference in unit tests; **not** Enzyme-reverse
compatible and so is not on the NUTS gradient path.
"""
struct HCubature <: ClusterAlg
    reltol::Float64
    abstol::Float64
    z_bound::Float64
end
HCubature(; reltol::Real = 1e-9, abstol::Real = 1e-12, z_bound::Real = 8.0) =
    HCubature(reltol, abstol, z_bound)

# Gauss-Hermite nodes/weights cache. We pay the FastGaussQuadrature
# construction cost once per (alg.n) value; the cache lets the model use
# the default n = 40 with no per-call overhead.
const _GH_CACHE = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
function _gh_nodes(n::Int)
    haskey(_GH_CACHE, n) && return _GH_CACHE[n]
    nodes, weights = gausshermite(n)
    out = (sqrt(2) .* nodes, weights ./ sqrt(π))
    _GH_CACHE[n] = out
    return out
end

# LogNormal CDF, written explicitly in terms of erf to keep the integrand
# allocation-free and to remove a Distributions.jl struct construction
# from the Enzyme tape.
_lognormal_cdf(s, μ, σ) =
    s > 0 ? (one(s) + erf((log(s) - μ) / (σ * sqrt(2)))) / 2 : zero(s)

"""
    F_cluster(t, μ_inc, σ_inc, μ_δ, σ_δ; alg = GaussHermite())

Probability that the full source-to-secondary chain
`Inc(src) + δ + Inc(sec)` is no greater than `t`, with
`Inc(src), Inc(sec) ~ LogNormal(μ_inc, σ_inc)` (i.i.d.) and
`δ ~ Normal(μ_δ, σ_δ)`.

`alg` selects the quadrature back-end (see [`GaussHermite`](@ref) and
[`HCubature`](@ref)). The default is the Enzyme-reverse-safe
Gauss-Hermite path.
"""
F_cluster(t::Real, μ_inc::Real, σ_inc::Real, μ_δ::Real, σ_δ::Real;
          alg::ClusterAlg = GaussHermite()) =
    _f_cluster(t, μ_inc, σ_inc, μ_δ, σ_δ, alg)

function _f_cluster(t, μ_inc, σ_inc, μ_δ, σ_δ, alg::GaussHermite)
    T = float(promote_type(typeof(t), typeof(μ_inc), typeof(σ_inc),
                           typeof(μ_δ), typeof(σ_δ)))
    t <= 0 && return zero(T)
    nodes, weights = _gh_nodes(alg.n)
    acc = zero(T)
    @inbounds for i in eachindex(nodes)
        z1 = nodes[i]
        w1 = weights[i]
        inc_src = exp(μ_inc + σ_inc * z1)
        for j in eachindex(nodes)
            z2 = nodes[j]
            w2 = weights[j]
            δ  = μ_δ + σ_δ * z2
            s  = t - inc_src - δ
            acc += w1 * w2 * _lognormal_cdf(s, μ_inc, σ_inc)
        end
    end
    return acc
end

function _f_cluster(t, μ_inc, σ_inc, μ_δ, σ_δ, alg::HCubature)
    t <= 0 && return zero(float(t))
    p = (; t, μ_inc, σ_inc, μ_δ, σ_δ)
    zb = alg.z_bound
    prob = IntegralProblem(_f_cluster_integrand_z,
                           ([-zb, -zb], [zb, zb]), p)
    sol = solve(prob, HCubatureJL();
                reltol = alg.reltol, abstol = alg.abstol)
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

# Vector-input wrapper for AD testing via DifferentiationInterface.jl.
# `θ = [t, μ_inc, σ_inc, μ_δ, σ_δ]`.
"""
    F_cluster_vec(θ; alg = GaussHermite())

Vector-input wrapper of [`F_cluster`](@ref) for AD testing with
`DifferentiationInterface.jl`. `θ = [t, μ_inc, σ_inc, μ_δ, σ_δ]`.
"""
F_cluster_vec(θ; alg::ClusterAlg = GaussHermite()) =
    F_cluster(θ[1], θ[2], θ[3], θ[4], θ[5]; alg = alg)
