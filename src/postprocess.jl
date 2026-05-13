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

Compute summary MCMC diagnostics from a sampled chain.

Returns the largest split-rhat and smallest bulk ESS over all scalar
parameter entries together with the total number of NUTS divergences
(`numerical_error`).

# Arguments
- `chn`: a chain returned by `Turing.sample`.

# Returns
A named tuple `(rhat, ess, ndiv)` with the maximum rhat, minimum bulk
ESS, and divergence count.
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

Extract pooled posterior draws for each entry of a vector-valued parameter.

For a chain holding a vector parameter such as `T_inf`, `T_onset`, or
`log_R`, returns one pooled-across-chains sample vector per entry.

# Arguments
- `chn`: a chain returned by `Turing.sample`.
- `name`: the parameter name as a `Symbol`, e.g. `:log_R`.

# Returns
A `Vector{Vector{Float64}}` of length equal to the parameter dimension,
each element containing the pooled posterior draws for that entry.
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

Print and return a posterior summary of the fitted model.

Computes the pooled posterior of the population parameters and the derived
generation-interval / serial-interval marginals and pre-symptomatic
transmission probabilities, prints [`summary_table`](@ref) to `stdout`,
and returns the underlying samples that [`save_posterior`](@ref) and the
plotting functions in `plots.jl` consume.

# Arguments
- `chn`: a chain returned by `Turing.sample`, e.g. from [`joint_model`](@ref).

# Returns
A named tuple `(μ_inc, σ_inc, μ_δ, σ_δ, k, log_R_chain, mean_gi_si,
sd_gi_si, p_pre)` of pooled posterior samples and derived quantities.
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

Write the posterior samples to a CSV file.

Writes one row per posterior draw with columns for the population
parameters, the derived GI / SI mean and SD, the pre-symptomatic
transmission probabilities, and one `log_R_b` column per weekly knot.
Creates the parent directory if needed.

# Arguments
- `post`: the summary tuple returned by [`summarise`](@ref).
- `path`: output CSV path.

# Returns
The result of `CSV.write` (the path written), as a side effect of writing
the file.
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
