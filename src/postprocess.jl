## Posterior summaries, CSV output, and R(t) figure.

# Analytic mean and variance of SkewNormal(ξ, ω, α) so we don't pull in
# Statistics.mean / var (collides with Statistics.quantile we already use).
function SkewNormal_meanvar(ξ, ω, α)
    d = α / sqrt(1 + α^2)
    m = ξ + ω * d * sqrt(2 / π)
    v = ω^2 * (1 - 2 * d^2 / π)
    return m, v
end

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
    α_δ   = vec(collect(chn[:α_δ]))
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

    println("Transmission timing δ (SkewNormal; days from source onset to secondary infection)")
    _print_qci("μ_δ (location, d)", μ_δ)
    _print_qci("σ_δ (scale, d)",    σ_δ)
    _print_qci("α_δ (skewness)",    α_δ)
    # SkewNormal has no analytic CDF in Distributions.jl; use MC samples
    p_pre = Dict{Float64,Vector{Float64}}()
    nmc = 5000
    rng_pp = Random.MersenneTwister(0)
    sn_samples = [rand(rng_pp, SkewNormal(μ_δ[i], σ_δ[i], α_δ[i]), nmc) for i in eachindex(μ_δ)]
    for τ in (0.0, -1.0, -2.0)
        p_pre[τ] = [sum(s .< τ) / length(s) for s in sn_samples]
        _print_qci(@sprintf("P(δ < %.0f)", τ), p_pre[τ]; fmt = "%.3f")
    end
    println()

    # Generation interval and serial interval as derived population marginals
    # (GI = δ + Inc_source, SI = δ + Inc_secondary). δ ~ SkewNormal.
    δ_dist_mean = [SkewNormal_meanvar(μ_δ[i], σ_δ[i], α_δ[i])[1] for i in eachindex(μ_δ)]
    δ_dist_var  = [SkewNormal_meanvar(μ_δ[i], σ_δ[i], α_δ[i])[2] for i in eachindex(μ_δ)]
    mean_gi_si  = δ_dist_mean .+ mean_inc
    sd_gi_si    = sqrt.(δ_dist_var .+ var_inc)
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

    return (; μ_inc, σ_inc, μ_δ, σ_δ, α_δ, k = k_s, log_R_chain,
            mean_gi_si, sd_gi_si, p_pre)
end

# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# R(t) figure
# ---------------------------------------------------------------------------

# Spaghetti plot of R(t) over weekly bins: each thinned posterior draw is a
# horizontal line for each bin, with the line broken between bins so no
# vertical step connectors are drawn. Bin edges come from `BIN_EDGES`
# (data.jl); the first and last bins extend one bin-width past the listed
# edges. Saved as a PNG.
function plot_rt(post, path; n_draws_plot = 100, ymax = 4.0)
    log_R = post.log_R_chain
    n_draws = length(log_R[1])
    step    = max(1, n_draws ÷ n_draws_plot)
    idx     = 1:step:n_draws

    bin_width  = BIN_EDGES[2] - BIN_EDGES[1]
    left_edge  = vcat(BIN_EDGES[1] - bin_width, BIN_EDGES)
    right_edge = vcat(BIN_EDGES, BIN_EDGES[end] + bin_width)
    # Two points per bin (left, right) then a NaN-on-y placeholder so the
    # line is broken between bins (no vertical step connectors).
    xs = Date[]
    for b in eachindex(log_R)
        push!(xs, left_edge[b]); push!(xs, right_edge[b]); push!(xs, right_edge[b])
    end

    plt = plot(; ylims = (0.0, ymax),
                 xlabel = "Date", ylabel = "R(t)",
                 legend = false,
                 title  = "Time-varying reproduction number (weekly bins)")
    for d in idx
        ys = Float64[]
        for b in eachindex(log_R)
            r = exp(log_R[b][d])
            push!(ys, r); push!(ys, r); push!(ys, NaN)
        end
        plot!(plt, xs, ys; linecolor = :steelblue, linewidth = 1.6, alpha = 0.25)
    end
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
    α_δ     = vec(collect(chn[:α_δ]))

    medians = Float64[]
    for i in 1:data.N
        src = data.source_idx[i]
        src == 0 && continue
        push!(medians, quantile(t_inf[i] .- t_onset[src], 0.5))
    end

    μ_med = quantile(μ_δ, 0.5)
    σ_med = quantile(σ_δ, 0.5)
    α_med = quantile(α_δ, 0.5)

    plt = histogram(medians; bins = 15, normalize = :pdf,
                    label  = "per-pair posterior medians (N = $(length(medians)))",
                    xlabel = "δ (days from source onset)", ylabel = "density",
                    title  = "Per-pair δ vs fitted population SkewNormal")
    xs = range(μ_med - 4σ_med, μ_med + 4σ_med; length = 200)
    plot!(plt, xs, pdf.(SkewNormal(μ_med, σ_med, α_med), xs);
          linewidth = 2, label = "SkewNormal(μ_δ, σ_δ, α_δ) fitted")
    vline!(plt, [0.0]; linestyle = :dash, color = :grey, label = "source onset")
    mkpath(dirname(path))
    savefig(plt, path)
    return path
end

# Pairplot of the population parameters: histograms on the diagonal, scatter
# in the lower triangle, blank in the upper. Subsampled for speed.
function plot_pairplot(post, path; n_draws_plot = 1000)
    params = [
        ("mu_inc", post.μ_inc),
        ("sig_inc", post.σ_inc),
        ("mu_delta", post.μ_δ),
        ("sig_delta", post.σ_δ),
        ("k", post.k),
    ]
    n = length(params)
    step = max(1, length(params[1][2]) ÷ n_draws_plot)
    idx  = 1:step:length(params[1][2])

    panels = []
    for i in 1:n, j in 1:n
        if i == j
            p = histogram(params[i][2]; bins = 30, normalize = :pdf,
                          legend = false, framestyle = :semi,
                          color = :steelblue, linecolor = :steelblue,
                          xlabel = (i == n) ? params[j][1] : "",
                          ylabel = (j == 1) ? params[i][1] : "",
                          xguidefontsize = 8, yguidefontsize = 8,
                          xtickfontsize = 6, ytickfontsize = 6)
        elseif i > j
            p = scatter(params[j][2][idx], params[i][2][idx];
                        markersize = 1.0, markeralpha = 0.10,
                        markerstrokewidth = 0, color = :steelblue,
                        legend = false, framestyle = :semi,
                        xlabel = (i == n) ? params[j][1] : "",
                        ylabel = (j == 1) ? params[i][1] : "",
                        xguidefontsize = 8, yguidefontsize = 8,
                        xtickfontsize = 6, ytickfontsize = 6)
        else
            p = plot(; framestyle = :none, legend = false)
        end
        push!(panels, p)
    end
    plt = plot(panels...; layout = (n, n), size = (1000, 1000),
                          plot_title = "Posterior pairplot",
                          plot_titlefontsize = 12)
    mkpath(dirname(path))
    savefig(plt, path)
    return path
end

# Prior predictives. Sample (μ_inc, σ_inc, μ_δ, σ_δ) from their independent
# priors, then for each draw sample Inc, δ, GI = δ + Inc. Show the implied
# distributions before seeing any data.
function plot_prior_predictives(path; n = 5000, rng = Random.MersenneTwister(0))
    μ_inc = rand(rng, Normal(3.0, 0.5), n)
    σ_inc = abs.(rand(rng, Normal(0.0, 0.5), n))
    μ_δ   = rand(rng, Normal(0.0, 5.0), n)
    σ_δ   = abs.(rand(rng, Normal(0.0, 1.0), n))
    inc_s = [rand(rng, LogNormal(μ_inc[i], σ_inc[i])) for i in 1:n]
    δ_s   = [rand(rng, Normal(μ_δ[i], σ_δ[i]))       for i in 1:n]
    gi_s  = δ_s .+ inc_s

    p_inc = histogram(inc_s; bins = 100, normalize = :pdf,
                      title = "Inc (prior)", xlabel = "days",
                      xlims = (0, 80), legend = false, color = :steelblue)
    p_del = histogram(δ_s; bins = 100, normalize = :pdf,
                      title = "δ (prior)", xlabel = "days from source onset",
                      xlims = (-25, 25), legend = false, color = :steelblue)
    p_gi  = histogram(gi_s; bins = 100, normalize = :pdf,
                      title = "GI / SI (prior)", xlabel = "days",
                      xlims = (-30, 80), legend = false, color = :steelblue)
    plt = plot(p_inc, p_del, p_gi; layout = (1, 3), size = (1500, 400))
    mkpath(dirname(path))
    savefig(plt, path)
    return path
end

# Posterior predictive: sample new pairs from the fitted (μ_*, σ_*) and
# overlay with the per-pair posterior medians from the observed pairs.
# If predicted and observed match in shape, the population summary is
# consistent with the data; deviations point at residual structure.
function plot_posterior_predictions(chn, data, path;
                                    n_per_draw = 50,
                                    rng = Random.MersenneTwister(123))
    μ_inc = vec(collect(chn[:μ_inc]))
    σ_inc = vec(collect(chn[:σ_inc]))
    μ_δ   = vec(collect(chn[:μ_δ]))
    σ_δ   = vec(collect(chn[:σ_δ]))
    α_δ   = vec(collect(chn[:α_δ]))
    n_draws = length(μ_inc)

    δ_pred  = Float64[]
    si_pred = Float64[]
    for i in 1:n_draws
        for _ in 1:n_per_draw
            d = rand(rng, SkewNormal(μ_δ[i], σ_δ[i], α_δ[i]))
            c = rand(rng, LogNormal(μ_inc[i], σ_inc[i]))
            push!(δ_pred,  d)
            push!(si_pred, d + c)
        end
    end

    t_inf   = vector_chain(chn, :T_inf)
    t_onset = vector_chain(chn, :T_onset)
    obs_δ  = Float64[]
    obs_si = Float64[]
    for i in 1:data.N
        src = data.source_idx[i]
        src == 0 && continue
        push!(obs_δ,  quantile(t_inf[i]   .- t_onset[src], 0.5))
        push!(obs_si, quantile(t_onset[i] .- t_onset[src], 0.5))
    end

    p_δ = histogram(δ_pred; bins = 80, normalize = :pdf,
                    alpha = 0.5, color = :steelblue, label = "posterior predictive",
                    title = "δ: predictive vs observed", xlabel = "days")
    histogram!(p_δ, obs_δ; bins = 15, normalize = :pdf,
               alpha = 0.5, color = :darkorange,
               label = "observed per-pair medians (N = $(length(obs_δ)))")

    p_si = histogram(si_pred; bins = 80, normalize = :pdf,
                     alpha = 0.5, color = :steelblue, label = "posterior predictive",
                     title = "SI: predictive vs observed", xlabel = "days")
    histogram!(p_si, obs_si; bins = 15, normalize = :pdf,
               alpha = 0.5, color = :darkorange,
               label = "observed per-pair medians (N = $(length(obs_si)))")

    plt = plot(p_δ, p_si; layout = (1, 2), size = (1400, 450))
    mkpath(dirname(path))
    savefig(plt, path)
    return path
end

function save_posterior(post, path)
    df = DataFrame(μ_inc = post.μ_inc, σ_inc = post.σ_inc,
                   μ_δ = post.μ_δ,     σ_δ = post.σ_δ, α_δ = post.α_δ,
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
