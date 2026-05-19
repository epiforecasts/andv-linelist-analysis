## Posterior summaries and CSV output. Plotting lives in `plots.jl`.

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

# Flat vector of scalar stat values, one per (scalar) parameter entry.
function _scalar_stats(summary)
    out = Float64[]
    for p in FlexiChains.parameters(summary)
        v = summary[p]
        if v isa Number
            ismissing(v) && continue
            push!(out, Float64(v))
        else
            for x in skipmissing(vec(collect(v)))
                push!(out, Float64(x))
            end
        end
    end
    return out
end

function _num_divergences(chn)
    for e in FlexiChains.extras(chn)
        e.name === :numerical_error || continue
        return Int(sum(skipmissing(vec(chn[e]))))
    end
    return 0
end

"""
$(TYPEDSIGNATURES)

Return convergence diagnostics for `chn`: `(; rhat, ess, ndiv)` — the
maximum `R̂` across scalar parameter entries, the minimum bulk ESS,
and the divergent transition count.

# Arguments
- `chn`: FlexiChain returned by [`sample_fit`](@ref).
"""
function diagnostics(chn)
    rhats = _scalar_stats(FlexiChains.rhat(chn))
    esses = _scalar_stats(FlexiChains.ess(chn; kind = :bulk))
    return (; rhat = maximum(rhats), ess = minimum(esses),
        ndiv = _num_divergences(chn))
end

# ---------------------------------------------------------------------------
# Extracting vector parameters
# ---------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Return a vector of pooled posterior samples for each entry of a
vector-valued parameter (e.g. `:T_inf`, `:log_R`).

# Arguments
- `chn`: FlexiChain returned by [`sample_fit`](@ref).
- `name`: name of a vector-valued parameter in `chn`.
"""
function vector_chain(chn, name::Symbol)
    arr = chn[name, stack = true]
    N = size(arr, 3)
    return [vec(collect(arr[:, :, i])) for i in 1:N]
end

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

function qci(x; q = (0.025, 0.5, 0.975))
    (quantile(x, q[1]), quantile(x, q[2]), quantile(x, q[3]))
end

function _print_qci(label, x; fmt = "%.2f")
    lo, med, hi = qci(x)
    @eval @printf $("  %-30s " * fmt * " (95%% CrI " * fmt * " – " * fmt *
                    ")\n") $label $med $lo $hi
end

"""
$(TYPEDSIGNATURES)

Build the named-tuple of posterior draws consumed by [`save_posterior`](@ref).
Call [`summary_table`](@ref) directly to print the headline summary.

# Arguments
- `chn`: FlexiChain returned by [`sample_fit`](@ref).
"""
function summarise(chn)
    μ_inc = vec(collect(chn[:μ_inc]));
    σ_inc = vec(collect(chn[:σ_inc]))
    μ_δ = vec(collect(chn[:μ_δ]));
    σ_δ = vec(collect(chn[:σ_δ]))
    k_s = vec(collect(chn[:k]))

    mean_inc = exp.(μ_inc .+ σ_inc .^ 2 ./ 2)
    var_inc = exp.(2μ_inc .+ σ_inc .^ 2) .* (exp.(σ_inc .^ 2) .- 1)
    mean_gi_si = μ_δ .+ mean_inc
    sd_gi_si = sqrt.(σ_δ .^ 2 .+ var_inc)

    p_pre = Dict{Float64, Vector{Float64}}()
    for τ in (0.0, -1.0, -2.0)
        p_pre[τ] = [cdf(Normal(μ_δ[i], σ_δ[i]), τ) for i in eachindex(μ_δ)]
    end

    log_R_chain = vector_chain(chn, :log_R)

    return (; μ_inc, σ_inc, μ_δ, σ_δ, k = k_s, log_R_chain,
        mean_gi_si, sd_gi_si, p_pre)
end

"""
$(TYPEDSIGNATURES)

Quantile band for R(t) across knots from one posterior.
Returns a `DataFrame` with one row per knot and columns `bin`, `lo`, `med`,
`hi` taken from the three quantiles `q` of `exp.(post.log_R_chain[b])`.
When `t0` is supplied, an additional `:date` column maps each knot to its
calendar date via `t0 + Day(edges[b])` (using `prepare_rt_edges(t0)` when
`edges` is not given).

# Arguments
- `post`: posterior summary named tuple from [`summarise`](@ref).

# Keyword Arguments
- `q`: three-tuple of lower, central, and upper quantiles.
- `t0`: time origin (`Date`) for the fit. When supplied, the returned
  frame includes a `:date` column for plotting on a calendar-date axis.
- `edges`: knot positions in days from `t0`, as returned by
  [`prepare_rt_edges`](@ref). Defaults to `prepare_rt_edges(t0)`.
"""
function rt_band(post; q = (0.1, 0.5, 0.9),
        t0::Union{Nothing, Date} = nothing,
        edges::Union{Nothing, AbstractVector} = nothing)
    chain = post.log_R_chain
    bins = collect(eachindex(chain))
    lo = [quantile(exp.(chain[b]), q[1]) for b in bins]
    med = [quantile(exp.(chain[b]), q[2]) for b in bins]
    hi = [quantile(exp.(chain[b]), q[3]) for b in bins]
    df = DataFrame(bin = bins, lo = lo, med = med, hi = hi)
    if t0 !== nothing
        e = edges === nothing ? prepare_rt_edges(t0) : edges
        length(e) >= length(bins) ||
            error("length(edges)=$(length(e)) < n_knots=$(length(bins))")
        df.date = Date[t0 + Day(round(Int, e[b])) for b in bins]
    end
    return df
end

"""
$(TYPEDSIGNATURES)

Extract simulated offspring counts `Zobs` from the NamedTuple returned by
`rand(rng, fix(joint_model(...), truth))`. `case_model` declares
`Z[i] ~ safe_nb(...)`, so each sampled `Z[i]` appears under a VarName
whose string form is `"Z[i]"`. Matching by string form keeps the helper
decoupled from AbstractPPL's internal VarName API.

Used by simulation-based recovery checks (see the sim-recovery
walkthrough and `test/test_joint_recovery.jl`) to round-trip a fixed
truth through the joint generative model and back into a refit.

# Arguments
- `sim`: NamedTuple of sampled values returned by `rand` on a `fix`-ed
  `joint_model`.
- `N`: number of cases in the line list (`d.N`); used to allocate the
  result and to validate that every `Z[i]` was sampled.
"""
function extract_simulated_Zobs(sim, N::Integer)
    Z = Vector{Int}(undef, N)
    seen = falses(N)
    for (k, v) in pairs(sim)
        ks = string(k)
        if startswith(ks, "Z[") && endswith(ks, "]")
            idx = parse(Int, ks[3:(end - 1)])
            Z[idx] = Int(v)
            seen[idx] = true
        end
    end
    all(seen) ||
        error("rand() returned only $(count(seen))/$N Z values; " *
              "indices missing: $(findall(!, seen))")
    return Z
end

"""
$(TYPEDSIGNATURES)

Single-row summary of predictive `samples`.
Returns a `NamedTuple` `(med, lo, hi, mean)` with `med`, `lo`, `hi`
taken from the quantiles `q` and rounded to `Int`, and `mean` the raw
sample mean.

# Arguments
- `samples`: vector of predictive draws (typically integer counts).

# Keyword Arguments
- `q`: three-tuple of lower, central, and upper quantiles.
"""
function summarise_predictive(samples::AbstractVector{<:Real};
        q = (0.1, 0.5, 0.9))
    med = Int(round(quantile(samples, q[2])))
    lo = Int(round(quantile(samples, q[1])))
    hi = Int(round(quantile(samples, q[3])))
    return (; med, lo, hi, mean = mean(samples))
end

# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Write the posterior summary `post` (as returned by [`summarise`](@ref))
to a CSV at `path`, one column per scalar parameter plus one column per
`log_R` knot.

# Arguments
- `post`: posterior named tuple from [`summarise`](@ref).
- `path`: output CSV path.
"""
function save_posterior(post, path)
    df = DataFrame(μ_inc = post.μ_inc, σ_inc = post.σ_inc,
        μ_δ = post.μ_δ, σ_δ = post.σ_δ,
        k = post.k,
        mean_gi_si = post.mean_gi_si, sd_gi_si = post.sd_gi_si,
        p_pre_0 = post.p_pre[0.0],
        p_pre_1 = post.p_pre[-1.0],
        p_pre_2 = post.p_pre[-2.0])
    for b in eachindex(post.log_R_chain)
        df[!, Symbol("log_R_$b")] = post.log_R_chain[b]
    end
    mkpath(dirname(path))
    CSV.write(path, df)
end
