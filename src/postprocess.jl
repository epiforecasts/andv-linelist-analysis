## Posterior summaries, CSV output, and R(t) figure.

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

function diagnostics(chn)
    rhats = _scalar_stats(MCMCChains.rhat(chn))
    esses = _scalar_stats(MCMCChains.ess(chn; kind = :bulk))
    return (; rhat = maximum(rhats), ess = minimum(esses), ndiv = _num_divergences(chn))
end

# ---------------------------------------------------------------------------
# Extracting vector parameters
# ---------------------------------------------------------------------------

# Returns a vector of pooled samples for each entry of a vector-valued
# parameter (T_inf, offset, log_R, ...).
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

function summarise(chn)
    d = diagnostics(chn)
    @printf "Diagnostics: max rhat = %.3f, min ess_bulk = %.0f, divergences = %d\n\n" d.rhat d.ess d.ndiv

    μ_inc = vec(collect(chn[:μ_inc])); σ_inc = vec(collect(chn[:σ_inc]))
    μ_δ   = vec(collect(chn[:μ_δ]));   σ_δ   = vec(collect(chn[:σ_δ]))
    k_s   = vec(collect(chn[:k]))

    mean_inc = exp.(μ_inc .+ σ_inc .^ 2 ./ 2)
    var_inc  = exp.(2μ_inc .+ σ_inc .^ 2) .* (exp.(σ_inc .^ 2) .- 1)
    q95_inc  = [quantile(LogNormal(μ_inc[i], σ_inc[i]), 0.95) for i in eachindex(μ_inc)]
    q99_inc  = [quantile(LogNormal(μ_inc[i], σ_inc[i]), 0.99) for i in eachindex(μ_inc)]
    println("Incubation period")
    _print_qci("mean (d)", mean_inc)
    _print_qci("95th percentile (d)", q95_inc)
    _print_qci("99th percentile (d)", q99_inc)
    println()

    println("Transmission timing δ (days from source onset to secondary infection)")
    _print_qci("μ_δ (d)", μ_δ)
    _print_qci("σ_δ (d)", σ_δ)
    p_pre = Dict{Float64,Vector{Float64}}()
    for τ in (0.0, -1.0, -2.0)
        p_pre[τ] = [cdf(Normal(μ_δ[i], σ_δ[i]), τ) for i in eachindex(μ_δ)]
        _print_qci(@sprintf("P(δ < %.0f)", τ), p_pre[τ]; fmt = "%.3f")
    end
    println()

    # Generation interval and serial interval as derived population marginals
    # (GI = δ + Inc_source, SI = δ + Inc_secondary). They coincide in mean and
    # SD because Inc_source and Inc_secondary are exchangeable in this model.
    mean_gi_si = μ_δ .+ mean_inc
    sd_gi_si   = sqrt.(σ_δ .^ 2 .+ var_inc)
    println("Generation interval / serial interval (derived: δ + Inc)")
    _print_qci("mean (d)", mean_gi_si)
    _print_qci("SD (d)",   sd_gi_si)
    println()

    println("Offspring distribution")
    _print_qci("dispersion k", k_s)
    println()

    log_R_chain = vector_chain(chn, :log_R)
    labels = bin_labels()
    println("R(t) by bin")
    for b in eachindex(log_R_chain)
        _print_qci(labels[b], exp.(log_R_chain[b]))
    end

    return (; μ_inc, σ_inc, μ_δ, σ_δ, k = k_s, log_R_chain,
            mean_gi_si, sd_gi_si, p_pre)
end

# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# R(t) figure
# ---------------------------------------------------------------------------

# Step plot of the median posterior R(t) by bin with a 95% credible-interval
# ribbon. Bin edges come from `BIN_EDGES` (data.jl); the first and last bins
# extend one bin-width past the listed edges. Saved as a PNG.
function plot_rt(post, path)
    log_R = post.log_R_chain
    medians = [quantile(exp.(log_R[b]), 0.5)   for b in eachindex(log_R)]
    los     = [quantile(exp.(log_R[b]), 0.025) for b in eachindex(log_R)]
    his     = [quantile(exp.(log_R[b]), 0.975) for b in eachindex(log_R)]

    bin_width = BIN_EDGES[2] - BIN_EDGES[1]
    left_edge  = vcat(BIN_EDGES[1] - bin_width, BIN_EDGES)
    right_edge = vcat(BIN_EDGES, BIN_EDGES[end] + bin_width)
    # Build step coordinates: each bin contributes two points at its edges.
    xs = Date[]; ys_med = Float64[]; ys_lo = Float64[]; ys_hi = Float64[]
    for b in eachindex(medians)
        push!(xs, left_edge[b]);  push!(ys_med, medians[b]); push!(ys_lo, los[b]); push!(ys_hi, his[b])
        push!(xs, right_edge[b]); push!(ys_med, medians[b]); push!(ys_lo, los[b]); push!(ys_hi, his[b])
    end

    plt = plot(xs, ys_med;
               ribbon = (ys_med .- ys_lo, ys_hi .- ys_med),
               fillalpha = 0.2, linewidth = 2,
               ylims = (0.0, 5.0),
               xlabel = "Date", ylabel = "R(t)",
               legend = false,
               title  = "Time-varying reproduction number (weekly bins)")
    hline!(plt, [1.0]; linestyle = :dash, color = :grey)
    mkpath(dirname(path))
    savefig(plt, path)
    return path
end

# Sense-check: compare the per-pair posterior of δ to the fitted population
# Normal(μ_δ, σ_δ). For each sourced pair, take the posterior of
# δ_pair = T_inf[secondary] − T_onset[source] and reduce to its median;
# then plot the histogram of those 33 per-pair medians with the population
# density overlaid. If the histogram spread tracks σ_δ and the centre tracks
# μ_δ, the population-level summary is faithful to the per-pair posterior.
function plot_delta_sense_check(chn, data, path)
    t_inf   = vector_chain(chn, :T_inf)
    t_onset = vector_chain(chn, :T_onset)
    μ_δ     = vec(collect(chn[:μ_δ]))
    σ_δ     = vec(collect(chn[:σ_δ]))

    medians = Float64[]
    for i in 1:data.N
        src = data.source_idx[i]
        src == 0 && continue
        push!(medians, quantile(t_inf[i] .- t_onset[src], 0.5))
    end

    μ_med = quantile(μ_δ, 0.5)
    σ_med = quantile(σ_δ, 0.5)

    plt = histogram(medians; bins = 15, normalize = :pdf,
                    label  = "per-pair posterior medians (N = $(length(medians)))",
                    xlabel = "δ (days from source onset)", ylabel = "density",
                    title  = "Per-pair δ vs fitted population Normal")
    xs = range(μ_med - 4σ_med, μ_med + 4σ_med; length = 200)
    plot!(plt, xs, pdf.(Normal(μ_med, σ_med), xs);
          linewidth = 2, label = "Normal(μ_δ, σ_δ) fitted")
    vline!(plt, [0.0]; linestyle = :dash, color = :grey, label = "source onset")
    mkpath(dirname(path))
    savefig(plt, path)
    return path
end

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
