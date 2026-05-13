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
"""
function diagnostics(chn)
    rhats = _scalar_stats(FlexiChains.rhat(chn))
    esses = _scalar_stats(FlexiChains.ess(chn; kind = :bulk))
    return (; rhat = maximum(rhats), ess = minimum(esses), ndiv = _num_divergences(chn))
end

# ---------------------------------------------------------------------------
# Extracting vector parameters
# ---------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Return a vector of pooled posterior samples for each entry of a
vector-valued parameter (e.g. `:T_inf`, `:log_R`).
"""
function vector_chain(chn, name::Symbol)
    arr = chn[name, stack = true]
    N = size(arr, 3)
    return [vec(collect(arr[:, :, i])) for i in 1:N]
end

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

qci(x; q = (0.025, 0.5, 0.975)) = (quantile(x, q[1]), quantile(x, q[2]), quantile(x, q[3]))

function _print_qci(label, x; fmt = "%.2f")
    lo, med, hi = qci(x)
    @eval @printf $("  %-30s " * fmt * " (95%% CrI " * fmt * " – " * fmt * ")\n") $label $med $lo $hi
end

"""
$(TYPEDSIGNATURES)

Build the named-tuple of posterior draws consumed by [`save_posterior`](@ref)
and print the headline summary table via [`summary_table`](@ref).
"""
function summarise(chn)
    μ_inc = vec(collect(chn[:μ_inc])); σ_inc = vec(collect(chn[:σ_inc]))
    μ_δ   = vec(collect(chn[:μ_δ]));   σ_δ   = vec(collect(chn[:σ_δ]))
    k_s   = vec(collect(chn[:k]))

    mean_inc = exp.(μ_inc .+ σ_inc .^ 2 ./ 2)
    var_inc  = exp.(2μ_inc .+ σ_inc .^ 2) .* (exp.(σ_inc .^ 2) .- 1)
    mean_gi_si = μ_δ .+ mean_inc
    sd_gi_si   = sqrt.(σ_δ .^ 2 .+ var_inc)

    p_pre = Dict{Float64,Vector{Float64}}()
    for τ in (0.0, -1.0, -2.0)
        p_pre[τ] = [cdf(Normal(μ_δ[i], σ_δ[i]), τ) for i in eachindex(μ_δ)]
    end

    log_R_chain = vector_chain(chn, :log_R)

    show(stdout, "text/plain", summary_table(chn))
    println()

    return (; μ_inc, σ_inc, μ_δ, σ_δ, k = k_s, log_R_chain,
            mean_gi_si, sd_gi_si, p_pre)
end

# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Write the posterior summary `post` (as returned by [`summarise`](@ref))
to a CSV at `path`, one column per scalar parameter plus one column per
`log_R` knot.
"""
function save_posterior(post, path)
    df = DataFrame(μ_inc = post.μ_inc, σ_inc = post.σ_inc,
                   μ_δ = post.μ_δ,     σ_δ = post.σ_δ,
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
