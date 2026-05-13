## Package-level plotting and summary-table helpers shared by the analysis
## walkthrough and the CLI. Every function returns a Makie `Figure` so the
## caller decides whether to render inline or write to disk. Plotting uses
## Makie + AlgebraOfGraphics; a Makie backend (e.g. CairoMakie) must be
## loaded at the call site to render or save figures.

# Apply a consistent theme to every figure produced here without mutating
# the user's global Makie theme.
_default_theme() = merge(theme_latexfonts(), Theme(fontsize = 12))

_with_theme(f) = with_theme(f, _default_theme())

"""
    plot_data(ll)

Two-panel view of the raw line list: epicurve by ISO week of onset (left)
and exposure windows against onset dates (right). Returns a `Makie.Figure`.
"""
function plot_data(ll)
    weekly = @chain DataFrame(week = ll.onset_date) begin
        @rtransform(:week = :week - Day(dayofweek(:week) - 1))
        @by(:week, :cases = length(:week))
        @orderby(:week)
    end

    sourced = @chain DataFrame(
        exposure_lower = ll.exposure_lower,
        exposure_upper = ll.exposure_upper,
        onset_date     = ll.onset_date) begin
        @rsubset(!ismissing(:exposure_lower))
        @orderby(:onset_date)
    end
    sourced[!, :idx] = 1:nrow(sourced)

    return _with_theme() do
        fig = Figure(; size = (1200, 450))

        # Epicurve via AlgebraOfGraphics: bar of weekly cases.
        weekly_aog = DataFrame(
            week_int = Dates.value.(weekly.week),
            cases    = weekly.cases,
        )
        epi_plot = data(weekly_aog) *
                   mapping(:week_int => "Week of onset",
                           :cases    => "Cases") *
                   visual(BarPlot, color = :steelblue)
        ag1 = draw!(fig[1, 1], epi_plot;
                    axis = (title = "Epicurve (weekly)",))
        ax1 = only(ag1).axis
        # Replace integer day-of-epoch ticks with the actual dates.
        let xs = Dates.value.(weekly.week)
            n = length(xs)
            keep = unique(round.(Int, range(1, n; length = min(n, 6))))
            ax1.xticks = (xs[keep], string.(weekly.week[keep]))
            ax1.xticklabelrotation = π / 6
        end

        ax2 = Axis(fig[1, 2]; xlabel = "Date",
                   ylabel = "Case (ordered by onset)",
                   title = "Exposure windows and onset")
        for r in eachrow(sourced)
            lines!(ax2,
                   [Dates.value(r.exposure_lower),
                    Dates.value(r.exposure_upper)],
                   [r.idx, r.idx];
                   color = :steelblue, linewidth = 3)
        end
        # Most exposure windows are a single day, so the line segment is
        # invisible; render a point at the window so 1-day cases show up.
        scatter!(ax2, Dates.value.(sourced.exposure_lower), sourced.idx;
                 color = :steelblue, markersize = 5,
                 label = "Exposure window / point")
        scatter!(ax2, Dates.value.(sourced.onset_date), sourced.idx;
                 color = :darkorange, markersize = 6, label = "Onset")
        axislegend(ax2; position = :rt, merge = true)
        let dts = sort(unique(vcat(sourced.exposure_lower,
                                   sourced.exposure_upper,
                                   sourced.onset_date)))
            n = length(dts)
            keep = unique(round.(Int, range(1, n; length = min(n, 6))))
            ax2.xticks = (Dates.value.(dts[keep]), string.(dts[keep]))
            ax2.xticklabelrotation = π / 6
        end

        fig
    end
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

# Add an inferred-density ribbon + predictive-sample histogram to `ax`. The
# histogram is normalised to a PDF so the two layers share a y-axis.
function _pp_panel!(ax, xs, dists, samples; bins = 50)
    med, lo, hi = _density_band(xs, dists)
    band!(ax, xs, lo, hi; color = (:steelblue, 0.25))
    inferred = lines!(ax, xs, med; color = :steelblue, linewidth = 2)
    pred = hist!(ax, samples; bins = bins, normalization = :pdf,
                 color = (:darkorange, 0.45),
                 strokecolor = :darkorange, strokewidth = 0.5)
    return inferred, pred
end

# Predictive-sample histogram only; used for GI / SI where no closed-form
# parametric density is available without further approximation.
function _pp_panel_hist_only!(ax, samples; bins = 50)
    return hist!(ax, samples; bins = bins, normalization = :pdf,
                 color = (:darkorange, 0.45),
                 strokecolor = :darkorange, strokewidth = 0.5)
end

"""
    plot_posterior_predictive(chn; rng = Random.MersenneTwister(1))

Two-by-two panel of posterior-predictive distributions for incubation
period, transmission timing δ, generation interval, and serial interval.
Inc and δ panels overlay the posterior over the parametric density
(median PDF with a 95% pointwise ribbon across draws) and a histogram of
one predictive sample per draw. GI and SI show the predictive-sample
histogram only. Returns a `Makie.Figure`.
"""
function plot_posterior_predictive(chn; rng = Random.MersenneTwister(1))
    samples = _ppc_frame(chn; rng)

    inc_dists = LogNormal.(samples.μ_inc, samples.σ_inc)
    δ_dists   = Normal.(samples.μ_δ, samples.σ_δ)

    return _with_theme() do
        fig = Figure(; size = (1000, 700))

        # GI/SI density not plotted: moment-matched Normal is a poor fit;
        # revisit with a proper KDE later.
        parametric_panels = [
            ("Incubation period",     range(0.5, 70.0; length = 200),
             inc_dists, vcat(samples.inc_src, samples.inc_sec), "days"),
            ("Transmission timing δ", range(-5.0, 5.0; length = 200),
             δ_dists, samples.δ, "days from source onset"),
        ]
        hist_only_panels = [
            ("Generation interval", samples.gi, "days"),
            ("Serial interval",     samples.si, "days"),
        ]

        local last_inferred = nothing
        local last_pred = nothing
        for (k, (title, xs, dists, samps, xlabel)) in enumerate(parametric_panels)
            row, col = fldmod1(k, 2)
            ax = Axis(fig[row, col]; title = title, xlabel = xlabel,
                      ylabel = "density",
                      titlesize = 11)
            last_inferred, last_pred = _pp_panel!(ax, xs, dists, samps)
        end
        for (k, (title, samps, xlabel)) in enumerate(hist_only_panels)
            row, col = fldmod1(k + length(parametric_panels), 2)
            ax = Axis(fig[row, col]; title = title, xlabel = xlabel,
                      ylabel = "density",
                      titlesize = 11)
            last_pred = _pp_panel_hist_only!(ax, samps)
        end

        Legend(fig[3, 1:2],
               [last_inferred, last_pred],
               ["inferred (median + 95%)", "predictive samples"];
               orientation = :horizontal, framevisible = false,
               tellheight = true, tellwidth = false)
        rowsize!(fig.layout, 3, Auto(0.1))
        fig
    end
end

"""
    plot_rt(chn; n_draws_plot = 100, ymax = 4.0)

Spaghetti plot of R(t) over weekly bins. Each thinned posterior draw is a
horizontal segment per bin with the line broken between bins so no vertical
step connectors are drawn. Bin edges come from `BIN_EDGES` (data.jl); the
first and last bins extend one bin-width past the listed edges. Returns a
`Makie.Figure`.
"""
function plot_rt(chn; n_draws_plot::Int = 100, ymax::Real = 4.0)
    log_R = vector_chain(chn, :log_R)
    n_draws = length(log_R[1])
    step    = max(1, n_draws ÷ n_draws_plot)
    idx     = 1:step:n_draws

    bin_width  = BIN_EDGES[2] - BIN_EDGES[1]
    left_edge  = vcat(BIN_EDGES[1] - bin_width, BIN_EDGES)
    right_edge = vcat(BIN_EDGES, BIN_EDGES[end] + bin_width)

    xs = Float64[]
    for b in eachindex(log_R)
        push!(xs, Dates.value(left_edge[b]))
        push!(xs, Dates.value(right_edge[b]))
        push!(xs, Dates.value(right_edge[b]))
    end

    return _with_theme() do
        fig = Figure(; size = (1000, 500))
        ax = Axis(fig[1, 1];
                  xlabel = "Date", ylabel = "R(t)",
                  title  = "Time-varying reproduction number (weekly bins)",
                  limits = (nothing, (0.0, ymax)))
        for d in idx
            ys = Float64[]
            for b in eachindex(log_R)
                r = exp(log_R[b][d])
                push!(ys, r); push!(ys, r); push!(ys, NaN)
            end
            lines!(ax, xs, ys;
                   color = (:steelblue, 0.25), linewidth = 1.6)
        end
        hlines!(ax, [1.0]; color = :grey, linestyle = :dash)

        # Date-formatted x ticks.
        edge_dates = sort(unique(vcat(left_edge, right_edge)))
        n = length(edge_dates)
        keep = unique(round.(Int, range(1, n; length = min(n, 7))))
        ax.xticks = (Dates.value.(edge_dates[keep]),
                     string.(edge_dates[keep]))
        ax.xticklabelrotation = π / 6
        fig
    end
end

"""
    plot_delta_sense_check(chn, data)

Sense-check the per-pair posterior of δ against the fitted population
`Normal(μ_δ, σ_δ)`. For each sourced pair, take the posterior of
`δ_pair = T_inf[secondary] − T_onset[source]` and reduce to its median; then
plot the histogram of those per-pair medians with the population density
overlaid. Returns a `Makie.Figure`.
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

    return _with_theme() do
        fig = Figure(; size = (900, 500))
        ax = Axis(fig[1, 1];
                  xlabel = "δ (days from source onset)",
                  ylabel = "density",
                  title  = "Per-pair δ vs fitted population Normal")
        h_per_pair = hist!(ax, medians;
                           bins = 15, normalization = :pdf,
                           color = (:steelblue, 0.55),
                           strokecolor = :steelblue, strokewidth = 0.5)
        xs = range(μ_med - 4σ_med, μ_med + 4σ_med; length = 200)
        l_fit = lines!(ax, xs, pdf.(Normal(μ_med, σ_med), xs);
                       color = :darkorange, linewidth = 2)
        v_zero = vlines!(ax, [0.0]; color = :grey, linestyle = :dash)
        Legend(fig[1, 2],
               [h_per_pair, l_fit, v_zero],
               ["per-pair posterior medians (N = $(length(medians)))",
                "Normal(μ_δ, σ_δ) fitted",
                "source onset"];
               framevisible = false, tellwidth = true)
        fig
    end
end

"""
    plot_prior_predictives(; n = 5000, rng = Random.MersenneTwister(0))

Prior-predictive panel: histograms of Inc, δ, and GI/SI drawn from the
package's independent priors on `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`. Returns a
`Makie.Figure`.
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

    df = vcat(
        DataFrame(value = inc_s, panel = "Inc (prior)",
                  xlabel = "days"),
        DataFrame(value = δ_s,   panel = "δ (prior)",
                  xlabel = "days from source onset"),
        DataFrame(value = gi_s,  panel = "GI / SI (prior)",
                  xlabel = "days"),
    )
    # Clip each panel to its viewing window so histograms aren't squashed
    # by long tails.
    windows = Dict("Inc (prior)" => (0.0, 80.0),
                   "δ (prior)"   => (-25.0, 25.0),
                   "GI / SI (prior)" => (-30.0, 80.0))
    df = @chain df begin
        @rsubset(windows[:panel][1] <= :value <= windows[:panel][2])
    end

    panels = [
        ("Inc (prior)",       inc_s, "days",                  (0.0, 80.0)),
        ("δ (prior)",         δ_s,   "days from source onset", (-25.0, 25.0)),
        ("GI / SI (prior)",   gi_s,  "days",                  (-30.0, 80.0)),
    ]

    return _with_theme() do
        fig = Figure(; size = (1500, 400))
        for (k, (title, samps, xlabel, xlim)) in enumerate(panels)
            ax = Axis(fig[1, k]; title = title, xlabel = xlabel,
                      ylabel = "density",
                      limits = (xlim, nothing))
            in_window = filter(x -> xlim[1] <= x <= xlim[2], samps)
            hist!(ax, in_window; bins = 100, normalization = :pdf,
                  color = :steelblue)
        end
        fig
    end
end
