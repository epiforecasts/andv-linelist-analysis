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
bulk ESS, divergence count, and wall-clock sampling time in seconds. The
runtime is read from FlexiChains' per-chain `sampling_time` metadata; under
`MCMCThreads` chains run in parallel so the wall clock is approximated by
the maximum over chains. Returns `missing` for the runtime if the chain
carries no timing metadata.
"""
function diagnostics_table(chn)
    d = diagnostics(chn)
    times = collect(skipmissing(FlexiChains.sampling_time(chn)))
    runtime = isempty(times) ? missing : maximum(times)
    return DataFrame(
        rhat_max          = d.rhat,
        ess_min           = d.ess,
        divergences       = d.ndiv,
        runtime_seconds   = runtime,
    )
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
    plot_predictive_distributions(chn; rng = Random.MersenneTwister(1))

Two-by-two panel of the implied population distributions under the
posterior for incubation period, transmission timing δ, generation
interval, and serial interval. Each panel shows draws from
`p(y_new | data) = ∫ p(y_new | θ) p(θ | data) dθ`, i.e. what a new case
or transmission pair would look like under the fitted parameters.

This is *not* a posterior-predictive check against observed data; for
that, see `plot_z_ppc`, `plot_delta_sense_check`, and
`plot_inc_sense_check`.

Inc and δ panels overlay the parametric density (median PDF with a 95%
pointwise ribbon across draws) and a histogram of one predictive sample
per draw. GI and SI show the predictive-sample histogram only. Returns a
`Makie.Figure`.
"""
function plot_predictive_distributions(chn; rng = Random.MersenneTwister(1))
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

Spaghetti plot of R(t) over the weekly knots. Each thinned posterior draw
is a piecewise-linear trajectory through `(knot_date[b], exp(log_R[b]))`.
Knot dates come from `BIN_EDGES` (data.jl). Returns a `Makie.Figure`.
"""
function plot_rt(chn; n_draws_plot::Int = 100, ymax::Real = 4.0)
    log_R = vector_chain(chn, :log_R)
    n_draws = length(log_R[1])
    step    = max(1, n_draws ÷ n_draws_plot)
    idx     = 1:step:n_draws

    knot_dates = BIN_EDGES
    xs = Float64[Dates.value(d) for d in knot_dates]

    return _with_theme() do
        fig = Figure(; size = (1000, 500))
        ax = Axis(fig[1, 1];
                  xlabel = "Date", ylabel = "R(t)",
                  title  = "Time-varying reproduction number (weekly knots)",
                  limits = (nothing, (0.0, ymax)))
        for d in idx
            ys = [exp(log_R[b][d]) for b in eachindex(log_R)]
            lines!(ax, xs, ys;
                   color = (:steelblue, 0.25), linewidth = 1.6)
        end
        hlines!(ax, [1.0]; color = :grey, linestyle = :dash)

        # Date-formatted x ticks.
        n = length(knot_dates)
        keep = unique(round.(Int, range(1, n; length = min(n, 7))))
        ax.xticks = (Dates.value.(knot_dates[keep]),
                     string.(knot_dates[keep]))
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
    plot_z_ppc(chn, data; rng = Random.MersenneTwister(1), edges = bin_edges_day(data.t0))

Posterior-predictive check for the observed offspring counts `Zobs`. For
each posterior draw and each case `i`, draws a replicated offspring count
`Z_rep[i] ~ NegativeBinomial(k, k/(k+R_i))`, where `R_i` is `exp(log_R)`
interpolated at the posterior median of `T_inf[i]` using the same machinery
as the model. Compares the count of cases at each `Z` value (0, 1, 2, …)
between the observed line list and the replicated distribution.

The left panel is a rootogram-style bar chart: bars are observed
frequencies; points + 95% pointwise CrI lines are replicated frequencies.
The right panel shows posterior-predictive distributions of three discrete
test statistics (`sum(Z)`, `max(Z)`, `count(Z == 0)`) with the observed
values marked.

Returns a `Makie.Figure`.
"""
function plot_z_ppc(chn, data;
                    rng = Random.MersenneTwister(1),
                    edges = bin_edges_day(data.t0))
    k_draws = _draws(chn, :k)
    log_R   = vector_chain(chn, :log_R)
    t_inf   = vector_chain(chn, :T_inf)
    n_draws = length(k_draws)
    N       = data.N

    # Use the posterior median of T_inf[i] for each case so the replicated
    # R_i tracks where the model places each case in time. This loses
    # within-case T_inf uncertainty but keeps the PPC interpretable.
    t_inf_med = [quantile(t_inf[i], 0.5) for i in 1:N]

    # Per-draw R_i: interpolate log_R at each case's median T_inf.
    R_per_draw = Matrix{Float64}(undef, n_draws, N)
    for d_idx in 1:n_draws
        logR_d = [log_R[b][d_idx] for b in eachindex(log_R)]
        for i in 1:N
            lr = log_R_at(t_inf_med[i], edges, logR_d)
            R_per_draw[d_idx, i] = exp(clamp(lr, -50.0, 50.0))
        end
    end

    # One Z_rep[i] per posterior draw.
    Z_rep = Matrix{Int}(undef, n_draws, N)
    for d_idx in 1:n_draws
        k_d = k_draws[d_idx]
        for i in 1:N
            R_i = R_per_draw[d_idx, i]
            p   = k_d / (k_d + R_i)
            Z_rep[d_idx, i] = rand(rng, NegativeBinomial(k_d, p))
        end
    end

    Zobs = data.Zobs
    zmax_obs = maximum(Zobs)
    # Show observed range plus a small margin; cap to keep the plot legible.
    zmax = min(max(zmax_obs + 2, 6), 20)
    z_values = 0:zmax

    obs_counts = [count(==(z), Zobs) for z in z_values]
    rep_counts = Matrix{Int}(undef, n_draws, length(z_values))
    for d_idx in 1:n_draws
        for (j, z) in enumerate(z_values)
            rep_counts[d_idx, j] = count(==(z), view(Z_rep, d_idx, :))
        end
    end
    rep_med = [quantile(view(rep_counts, :, j), 0.5)   for j in eachindex(z_values)]
    rep_lo  = [quantile(view(rep_counts, :, j), 0.025) for j in eachindex(z_values)]
    rep_hi  = [quantile(view(rep_counts, :, j), 0.975) for j in eachindex(z_values)]

    # Aggregate test statistics per draw.
    sum_rep   = [sum(view(Z_rep, d_idx, :))     for d_idx in 1:n_draws]
    max_rep   = [maximum(view(Z_rep, d_idx, :)) for d_idx in 1:n_draws]
    zeros_rep = [count(==(0), view(Z_rep, d_idx, :)) for d_idx in 1:n_draws]
    sum_obs   = sum(Zobs)
    max_obs   = maximum(Zobs)
    zeros_obs = count(==(0), Zobs)

    return _with_theme() do
        fig = Figure(; size = (1200, 500))

        ax1 = Axis(fig[1, 1];
                   xlabel = "Offspring count Z",
                   ylabel = "Number of cases",
                   title  = "Z by value: observed vs replicated",
                   xticks = collect(z_values))
        b_obs = barplot!(ax1, collect(z_values), Float64.(obs_counts);
                         color = (:steelblue, 0.55),
                         strokecolor = :steelblue, strokewidth = 0.5)
        # Error bars for the replicated 95% CrI plus a median marker.
        rangebars!(ax1, collect(z_values), rep_lo, rep_hi;
                   color = :darkorange, whiskerwidth = 8)
        s_rep = scatter!(ax1, collect(z_values), rep_med;
                         color = :darkorange, markersize = 8)

        ax2 = Axis(fig[1, 2];
                   xlabel = "Test statistic value",
                   ylabel = "density",
                   title  = "Discrete test statistics (replicated vs observed)")
        # Three overlaid histograms on the same x; each statistic gets its
        # own colour with the observed value as a dashed vline.
        h_sum = hist!(ax2, sum_rep; bins = 30, normalization = :pdf,
                      color = (:steelblue, 0.45))
        v_sum = vlines!(ax2, [Float64(sum_obs)];
                        color = :steelblue, linestyle = :dash, linewidth = 2)
        h_max = hist!(ax2, max_rep; bins = 30, normalization = :pdf,
                      color = (:darkorange, 0.45))
        v_max = vlines!(ax2, [Float64(max_obs)];
                        color = :darkorange, linestyle = :dash, linewidth = 2)
        h_zero = hist!(ax2, zeros_rep; bins = 30, normalization = :pdf,
                       color = (:seagreen, 0.45))
        v_zero = vlines!(ax2, [Float64(zeros_obs)];
                         color = :seagreen, linestyle = :dash, linewidth = 2)

        Legend(fig[2, 1],
               [b_obs, s_rep],
               ["observed", "replicated (median + 95% CrI)"];
               orientation = :horizontal, framevisible = false,
               tellheight = true, tellwidth = false)
        Legend(fig[2, 2],
               [h_sum, h_max, h_zero],
               ["sum(Z)", "max(Z)", "count(Z = 0)"];
               orientation = :horizontal, framevisible = false,
               tellheight = true, tellwidth = false)
        rowsize!(fig.layout, 2, Auto(0.12))
        fig
    end
end

"""
    plot_inc_sense_check(chn, data; n_density_draws = 200)

Sense-check the per-case posterior of the incubation period against the
fitted population `LogNormal(μ_inc, σ_inc)`. For each case, takes the
posterior of `inc_i = T_onset[i] − T_inf[i]` and reduces to its median;
plots the histogram of those per-case medians with the median PDF (and
95% pointwise ribbon) of the population LogNormal overlaid. Returns a
`Makie.Figure`.
"""
function plot_inc_sense_check(chn, data; n_density_draws::Int = 200)
    t_inf   = vector_chain(chn, :T_inf)
    t_onset = vector_chain(chn, :T_onset)
    μ_inc   = _draws(chn, :μ_inc)
    σ_inc   = _draws(chn, :σ_inc)

    medians = [quantile(t_onset[i] .- t_inf[i], 0.5) for i in 1:data.N]

    # Thin posterior draws of (μ_inc, σ_inc) for the density ribbon.
    step = max(1, length(μ_inc) ÷ n_density_draws)
    idx  = 1:step:length(μ_inc)
    dists = [LogNormal(μ_inc[i], σ_inc[i]) for i in idx]

    upper = max(maximum(medians) * 1.2, 60.0)
    xs = range(0.5, upper; length = 200)
    med, lo, hi = _density_band(xs, dists)

    return _with_theme() do
        fig = Figure(; size = (900, 500))
        ax = Axis(fig[1, 1];
                  xlabel = "Incubation period (days)",
                  ylabel = "density",
                  title  = "Per-case Inc vs fitted population LogNormal")
        h_per_case = hist!(ax, medians;
                           bins = 15, normalization = :pdf,
                           color = (:steelblue, 0.55),
                           strokecolor = :steelblue, strokewidth = 0.5)
        b_fit = band!(ax, xs, lo, hi; color = (:darkorange, 0.25))
        l_fit = lines!(ax, xs, med; color = :darkorange, linewidth = 2)
        Legend(fig[1, 2],
               [h_per_case, l_fit, b_fit],
               ["per-case posterior medians (N = $(length(medians)))",
                "LogNormal(μ_inc, σ_inc) median PDF",
                "95% pointwise ribbon"];
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
