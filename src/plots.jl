## Package-level plotting and summary-table helpers shared by the analysis
## walkthrough and the CLI. Every function returns a figure object so the
## caller decides whether to render inline or write to disk.

"""
    plot_data(ll)

Two-panel view of the raw line list: epicurve by ISO week of onset (left)
and exposure windows against onset dates (right). Returns a Plots.jl `Plot`.
"""
function plot_data(ll)
    weekly = @chain DataFrame(week = ll.onset_date) begin
        @rtransform(:week = :week - Day(dayofweek(:week) - 1))
        @by(:week, :cases = length(:week))
        @orderby(:week)
    end

    epicurve = bar(weekly.week, weekly.cases;
                   legend = false, xlabel = "Week of onset",
                   ylabel = "Cases", title = "Epicurve (weekly)",
                   color = :steelblue, linecolor = :steelblue)

    sourced = @chain DataFrame(
        exposure_lower = ll.exposure_lower,
        exposure_upper = ll.exposure_upper,
        onset_date     = ll.onset_date) begin
        @rsubset(!ismissing(:exposure_lower))
        @orderby(:onset_date)
    end

    expo = plot(; xlabel = "Date", ylabel = "Case (ordered by onset)",
                  title = "Exposure windows and onset", legend = :topright)
    for (j, r) in enumerate(eachrow(sourced))
        plot!(expo, [r.exposure_lower, r.exposure_upper], [j, j];
              linewidth = 3, color = :steelblue, label = "")
        scatter!(expo, [r.onset_date], [j];
                 markersize = 3, color = :darkorange, label = "")
    end
    scatter!(expo, Date[], Int[]; color = :steelblue,
             markershape = :hline, label = "Exposure window")
    scatter!(expo, Date[], Int[]; color = :darkorange, label = "Onset")

    return plot(epicurve, expo; layout = (1, 2), size = (1200, 450))
end

# Per-draw scalar parameter as a Vector{Float64} via FlexiChains. Used to
# build summaries and predictive samples without going through MCMCChains.
_draws(chn, name::Symbol) = vec(collect(chn[name]))

"""
    summary_table(chn)

Posterior summary `DataFrame` for the headline quantities: incubation mean,
95th and 99th percentiles, transmission timing μ_δ / σ_δ, GI / SI mean and
SD, and Negative-Binomial dispersion k. Columns: `quantity`, `median`,
`lower_95`, `upper_95`.
"""
function summary_table(chn)
    μ_inc = _draws(chn, :μ_inc)
    σ_inc = _draws(chn, :σ_inc)
    μ_δ   = _draws(chn, :μ_δ)
    σ_δ   = _draws(chn, :σ_δ)
    k     = _draws(chn, :k)

    per_draw = DataFrame(; μ_inc, σ_inc, μ_δ, σ_δ, k)
    derived  = @chain per_draw begin
        @transform begin
            :mean_inc = exp.(:μ_inc .+ :σ_inc .^ 2 ./ 2)
            :var_inc  = (exp.(:σ_inc .^ 2) .- 1) .* exp.(2 .* :μ_inc .+ :σ_inc .^ 2)
        end
        @rtransform begin
            :q95_inc   = quantile(LogNormal(:μ_inc, :σ_inc), 0.95)
            :q99_inc   = quantile(LogNormal(:μ_inc, :σ_inc), 0.99)
            :gi_si_mean = :μ_δ + :mean_inc
            :gi_si_sd   = sqrt(:σ_δ ^ 2 + :var_inc)
        end
    end

    rows = [
        ("Incubation mean (d)",       derived.mean_inc),
        ("Incubation 95th pct (d)",   derived.q95_inc),
        ("Incubation 99th pct (d)",   derived.q99_inc),
        ("μ_δ (d from source onset)", derived.μ_δ),
        ("σ_δ (d)",                   derived.σ_δ),
        ("GI / SI mean (d)",          derived.gi_si_mean),
        ("GI / SI SD (d)",            derived.gi_si_sd),
        ("NB dispersion k",           derived.k),
    ]

    return DataFrame(
        quantity = first.(rows),
        median   = [quantile(x, 0.5)   for (_, x) in rows],
        lower_95 = [quantile(x, 0.025) for (_, x) in rows],
        upper_95 = [quantile(x, 0.975) for (_, x) in rows],
    )
end

"""
    diagnostics_table(chn)

Single-row `DataFrame` summarising sampler diagnostics: maximum R̂, minimum
bulk ESS, and divergence count.
"""
function diagnostics_table(chn)
    d = diagnostics(chn)
    return DataFrame(rhat_max = d.rhat, ess_min = d.ess, divergences = d.ndiv)
end

"""
    plot_pair(chn; thin = 2)

Corner plot of the population scalars `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`, `k`
via PairPlots.jl. Returns a Makie `Figure` (requires a Makie backend such
as CairoMakie loaded at the call site).
"""
function plot_pair(chn; thin::Int = 2)
    tbl = @chain DataFrame(
        μ_inc = _draws(chn, :μ_inc),
        σ_inc = _draws(chn, :σ_inc),
        μ_δ   = _draws(chn, :μ_δ),
        σ_δ   = _draws(chn, :σ_δ),
        k     = _draws(chn, :k)) begin
        @rtransform(:_keep = true)
        @select(:μ_inc, :σ_inc, :μ_δ, :σ_δ, :k)
    end
    return pairplot(tbl[1:thin:end, :])
end

# Build the inferred-density ribbon and posterior-predictive sample frames
# in one pass so both the median PDF and the sample histogram come from the
# same set of posterior draws.
function _ppc_frame(chn; rng = Random.MersenneTwister(1))
    μ_inc = _draws(chn, :μ_inc)
    σ_inc = _draws(chn, :σ_inc)
    μ_δ   = _draws(chn, :μ_δ)
    σ_δ   = _draws(chn, :σ_δ)
    n = length(μ_inc)

    # One source-incubation and one secondary-incubation per draw so GI and
    # SI inherit independent incubation realisations as in the data.
    return @chain DataFrame(; μ_inc, σ_inc, μ_δ, σ_δ) begin
        @rtransform begin
            :inc_src = rand(rng, LogNormal(:μ_inc, :σ_inc))
            :inc_sec = rand(rng, LogNormal(:μ_inc, :σ_inc))
            :δ       = rand(rng, Normal(:μ_δ, :σ_δ))
        end
        @rtransform begin
            :gi = :δ + :inc_src
            :si = :δ + :inc_sec
        end
    end
end

function _density_band(xs, dists)
    pdfs = [pdf.(d, xs) for d in dists]
    mat  = reduce(hcat, pdfs)
    med  = [quantile(view(mat, j, :), 0.5)   for j in eachindex(xs)]
    lo   = [quantile(view(mat, j, :), 0.025) for j in eachindex(xs)]
    hi   = [quantile(view(mat, j, :), 0.975) for j in eachindex(xs)]
    return med, lo, hi
end

function _pp_panel(title, xs, dists, samples, xlabel; bins = 50)
    med, lo, hi = _density_band(xs, dists)
    plt = plot(xs, med; ribbon = (med .- lo, hi .- med),
               color = :steelblue, fillalpha = 0.25, linewidth = 2,
               label = "inferred (median + 95%)",
               title = title, xlabel = xlabel, ylabel = "density",
               legend = :topright, titlefontsize = 9, legendfontsize = 6)
    histogram!(plt, samples; bins = bins, normalize = :pdf,
               color = :darkorange, linecolor = :darkorange,
               alpha = 0.45, label = "predictive samples")
    return plt
end

"""
    plot_posterior_predictive(chn; rng = Random.MersenneTwister(1))

Two-by-two panel of posterior-predictive distributions for incubation
period, transmission timing δ, generation interval, and serial interval.
Each panel overlays the posterior over the parametric density (median PDF
with a 95% pointwise ribbon across draws) and a histogram of one
predictive sample per draw. GI / SI use a Normal moment-match of
`Normal(μ_δ, σ_δ) + LogNormal(μ_inc, σ_inc)` for the inferred ribbon and
exact `δ + Inc` draws for the histogram. Returns a Plots.jl `Plot`.
"""
function plot_posterior_predictive(chn; rng = Random.MersenneTwister(1))
    samples = _ppc_frame(chn; rng)

    inc_dists = LogNormal.(samples.μ_inc, samples.σ_inc)
    δ_dists   = Normal.(samples.μ_δ, samples.σ_δ)

    # Moment-matched Normal for GI / SI inferred ribbon — keeps the doc
    # build cheap while still showing the posterior over the marginal.
    mean_inc = exp.(samples.μ_inc .+ samples.σ_inc .^ 2 ./ 2)
    var_inc  = (exp.(samples.σ_inc .^ 2) .- 1) .*
               exp.(2 .* samples.μ_inc .+ samples.σ_inc .^ 2)
    gisi_dists = Normal.(samples.μ_δ .+ mean_inc,
                         sqrt.(samples.σ_δ .^ 2 .+ var_inc))

    p_inc = _pp_panel("Incubation period",     range(0.5, 70.0; length = 200),
                      inc_dists, vcat(samples.inc_src, samples.inc_sec), "days")
    p_δ   = _pp_panel("Transmission timing δ", range(-5.0, 5.0; length = 200),
                      δ_dists, samples.δ, "days from source onset")
    p_gi  = _pp_panel("Generation interval",   range(0.5, 80.0; length = 200),
                      gisi_dists, samples.gi, "days")
    p_si  = _pp_panel("Serial interval",       range(0.5, 80.0; length = 200),
                      gisi_dists, samples.si, "days")

    return plot(p_inc, p_δ, p_gi, p_si; layout = (2, 2), size = (1000, 700))
end

"""
    plot_rt(chn; n_draws_plot = 100, ymax = 4.0)

Spaghetti plot of R(t) over weekly bins. Each thinned posterior draw is a
horizontal segment per bin with the line broken between bins so no vertical
step connectors are drawn. Bin edges come from `BIN_EDGES` (data.jl); the
first and last bins extend one bin-width past the listed edges. Returns a
Plots.jl `Plot`.
"""
function plot_rt(chn; n_draws_plot::Int = 100, ymax::Real = 4.0)
    log_R = vector_chain(chn, :log_R)
    n_draws = length(log_R[1])
    step    = max(1, n_draws ÷ n_draws_plot)
    idx     = 1:step:n_draws

    bin_width  = BIN_EDGES[2] - BIN_EDGES[1]
    left_edge  = vcat(BIN_EDGES[1] - bin_width, BIN_EDGES)
    right_edge = vcat(BIN_EDGES, BIN_EDGES[end] + bin_width)
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
    return plt
end

"""
    plot_delta_sense_check(chn, data)

Sense-check the per-pair posterior of δ against the fitted population
`Normal(μ_δ, σ_δ)`. For each sourced pair, take the posterior of
`δ_pair = T_inf[secondary] − T_onset[source]` and reduce to its median; then
plot the histogram of those per-pair medians with the population density
overlaid. Returns a Plots.jl `Plot`.
"""
function plot_delta_sense_check(chn, data)
    t_inf   = vector_chain(chn, :T_inf)
    t_onset = vector_chain(chn, :T_onset)
    μ_δ     = _draws(chn, :μ_δ)
    σ_δ     = _draws(chn, :σ_δ)

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
    return plt
end

"""
    plot_prior_predictives(; n = 5000, rng = Random.MersenneTwister(0))

Prior-predictive panel: histograms of Inc, δ, and GI/SI drawn from the
package's independent priors on `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`. Returns a
Plots.jl `Plot`.
"""
function plot_prior_predictives(; n::Int = 5000,
                                  rng = Random.MersenneTwister(0))
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
    return plot(p_inc, p_del, p_gi; layout = (1, 3), size = (1500, 400))
end
