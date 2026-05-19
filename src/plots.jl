## Package-level plotting and summary-table helpers shared by the analysis
## walkthrough and the CLI. Every function returns a Makie `Figure` so the
## caller decides whether to render inline or write to disk. CairoMakie is
## loaded as the default backend by the top-level module.

# Apply a consistent theme to every figure produced here without mutating
# the user's global Makie theme.
_default_theme() = merge(theme_latexfonts(), Theme(fontsize = 12))

_with_theme(f) = with_theme(f, _default_theme())

"""
$(TYPEDSIGNATURES)

Two-panel view of the raw line list: epicurve by ISO week of onset (left)
and exposure windows against onset dates (right). Returns a `Makie.Figure`.

# Arguments
- `ll`: line-list `DataFrame` with `onset_date`, `exposure_lower`, and
  `exposure_upper` columns.
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
        onset_date = ll.onset_date) begin
        @rsubset(!ismissing(:exposure_lower))
        @orderby(:onset_date)
    end
    sourced[!, :idx] = 1:nrow(sourced)

    return _with_theme() do
        fig = Figure(; size = (1200, 450))

        # Epicurve via AlgebraOfGraphics: bar of weekly cases.
        weekly_aog = DataFrame(
            week_int = Dates.value.(weekly.week),
            cases = weekly.cases
        )
        epi_plot = data(weekly_aog) *
                   mapping(:week_int => "Week of onset",
                       :cases => "Cases") *
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

        # Exposure-window panel via AoG. Each case becomes two long-form
        # rows (lower, upper) so a grouped Lines layer draws one horizontal
        # segment per case; scatter layers add the onset point and a marker
        # at the lower endpoint (most windows are a single day, so the line
        # is invisible).
        segs = @chain sourced begin
            @select(:idx,
                :lower=Dates.value.(:exposure_lower),
                :upper=Dates.value.(:exposure_upper))
        end
        segs_long = vcat(
            @select(segs, :idx, :date = :lower),
            @select(segs, :idx, :date = :upper)
        )
        scatter_df = @chain sourced begin
            @select(:idx,
                :exposure_lower=Dates.value.(:exposure_lower),
                :onset_date=Dates.value.(:onset_date))
        end

        exposure_seg = data(segs_long) *
                       mapping(:date => "Date",
                           :idx => "Case (ordered by onset)",
                           group = :idx => nonnumeric) *
                       visual(Lines, color = :steelblue, linewidth = 3)
        exposure_pt = data(scatter_df) *
                      mapping(:exposure_lower => "Date",
                          :idx => "Case (ordered by onset)") *
                      visual(Scatter, color = :steelblue, markersize = 5)
        onset_pt = data(scatter_df) *
                   mapping(:onset_date => "Date",
                       :idx => "Case (ordered by onset)") *
                   visual(Scatter, color = :darkorange, markersize = 6)

        ag2 = draw!(fig[1, 2], exposure_seg + exposure_pt + onset_pt;
            axis = (title = "Exposure windows and onset",))
        ax2 = only(ag2).axis
        # Hand-rolled legend: AoG won't auto-build one from constant-style
        # layers with no shared categorical mapping.
        legend_elems = [
            MarkerElement(color = :steelblue, marker = :circle,
                markersize = 8),
            MarkerElement(color = :darkorange, marker = :circle,
                markersize = 8)
        ]
        axislegend(ax2, legend_elems,
            ["Exposure window / point", "Onset"];
            position = :rt)
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
$(TYPEDSIGNATURES)

Posterior summary `DataFrame` for the headline quantities: incubation mean,
95th and 99th percentiles, transmission timing μ_δ / σ_δ, GI / SI mean and
SD, and Negative-Binomial dispersion k. Columns: `quantity`, `median`,
`lower_95`, `upper_95`.

# Arguments
- `chn`: a sampled chain from [`joint_model`](@ref).
"""
function summary_table(chn)
    μ_inc = _draws(chn, :μ_inc)
    σ_inc = _draws(chn, :σ_inc)
    μ_δ = _draws(chn, :μ_δ)
    σ_δ = _draws(chn, :σ_δ)
    k = _draws(chn, :k)

    per_draw = DataFrame(; μ_inc, σ_inc, μ_δ, σ_δ, k)
    derived = @chain per_draw begin
        @transform begin
            :mean_inc = exp.(:μ_inc .+ :σ_inc .^ 2 ./ 2)
            :var_inc = (exp.(:σ_inc .^ 2) .- 1) .*
                       exp.(2 .* :μ_inc .+ :σ_inc .^ 2)
        end
        @rtransform begin
            :q95_inc = quantile(LogNormal(:μ_inc, :σ_inc), 0.95)
            :q99_inc = quantile(LogNormal(:μ_inc, :σ_inc), 0.99)
            :gi_si_mean = :μ_δ + :mean_inc
            :gi_si_sd = sqrt(:σ_δ ^ 2 + :var_inc)
        end
    end

    rows = [
        ("Incubation mean (d)", derived.mean_inc),
        ("Incubation 95th pct (d)", derived.q95_inc),
        ("Incubation 99th pct (d)", derived.q99_inc),
        ("μ_δ (d from source onset)", derived.μ_δ),
        ("σ_δ (d)", derived.σ_δ),
        ("GI / SI mean (d)", derived.gi_si_mean),
        ("GI / SI SD (d)", derived.gi_si_sd),
        ("NB dispersion k", derived.k)
    ]

    return DataFrame(
        quantity = first.(rows),
        median = [quantile(x, 0.5) for (_, x) in rows],
        lower_95 = [quantile(x, 0.025) for (_, x) in rows],
        upper_95 = [quantile(x, 0.975) for (_, x) in rows]
    )
end

"""
$(TYPEDSIGNATURES)

Single-row `DataFrame` summarising sampler diagnostics: maximum R̂, minimum
bulk ESS, divergence count, and wall-clock sampling time in seconds. The
runtime is read from FlexiChains' per-chain `sampling_time` metadata; under
`MCMCThreads` chains run in parallel so the wall clock is approximated by
the maximum over chains. Returns `missing` for the runtime if the chain
carries no timing metadata.

# Arguments
- `chn`: a sampled chain from [`joint_model`](@ref).
"""
function diagnostics_table(chn)
    d = diagnostics(chn)
    times = collect(skipmissing(FlexiChains.sampling_time(chn)))
    runtime = isempty(times) ? missing : maximum(times)
    return DataFrame(
        rhat_max = d.rhat,
        ess_min = d.ess,
        divergences = d.ndiv,
        runtime_seconds = runtime
    )
end

"""
$(TYPEDSIGNATURES)

Corner plot of the population scalars `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`, `k`
via PairPlots.jl. Returns a Makie `Figure`.

# Arguments
- `chn`: a sampled chain from [`joint_model`](@ref).

# Keyword Arguments
- `thin`: stride applied to the posterior draws before plotting.
"""
function plot_pair(chn; thin::Int = 2)
    tbl = @chain DataFrame(
        μ_inc = _draws(chn, :μ_inc),
        σ_inc = _draws(chn, :σ_inc),
        μ_δ = _draws(chn, :μ_δ),
        σ_δ = _draws(chn, :σ_δ),
        k = _draws(chn, :k)) begin
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
    μ_δ = _draws(chn, :μ_δ)
    σ_δ = _draws(chn, :σ_δ)

    # One source-incubation and one secondary-incubation per draw so GI and
    # SI inherit independent incubation realisations as in the data.
    return @chain DataFrame(; μ_inc, σ_inc, μ_δ, σ_δ) begin
        @rtransform begin
            :inc_src = rand(rng, LogNormal(:μ_inc, :σ_inc))
            :inc_sec = rand(rng, LogNormal(:μ_inc, :σ_inc))
            :δ = rand(rng, Normal(:μ_δ, :σ_δ))
        end
        @rtransform begin
            :gi = :δ + :inc_src
            :si = :δ + :inc_sec
        end
    end
end

function _density_band(xs, dists)
    pdfs = [pdf.(d, xs) for d in dists]
    mat = reduce(hcat, pdfs)
    med = [quantile(view(mat, j, :), 0.5) for j in eachindex(xs)]
    lo = [quantile(view(mat, j, :), 0.025) for j in eachindex(xs)]
    hi = [quantile(view(mat, j, :), 0.975) for j in eachindex(xs)]
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
$(TYPEDSIGNATURES)

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

# Arguments
- `chn`: a sampled chain from [`joint_model`](@ref).

# Keyword Arguments
- `rng`: RNG used to draw one predictive sample per posterior draw.
"""
function plot_predictive_distributions(chn; rng = Random.MersenneTwister(1))
    samples = _ppc_frame(chn; rng)

    inc_dists = LogNormal.(samples.μ_inc, samples.σ_inc)
    δ_dists = Normal.(samples.μ_δ, samples.σ_δ)

    return _with_theme() do
        fig = Figure(; size = (1000, 700))

        # GI/SI density not plotted: moment-matched Normal is a poor fit;
        # revisit with a proper KDE later.
        parametric_panels = [
            ("Incubation period", range(0.5, 70.0; length = 200),
                inc_dists, vcat(samples.inc_src, samples.inc_sec), "days"),
            ("Transmission timing δ", range(-5.0, 5.0; length = 200),
                δ_dists, samples.δ, "days from source onset")
        ]
        hist_only_panels = [
            ("Generation interval", samples.gi, "days"),
            ("Serial interval", samples.si, "days")
        ]

        local last_inferred = nothing
        local last_pred = nothing
        for (k, (title, xs, dists, samps, xlabel)) in
            enumerate(parametric_panels)
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
$(TYPEDSIGNATURES)

Spaghetti plot of R(t) over the weekly knots. Each thinned posterior draw
is a piecewise-linear trajectory through `(knot_date[b], exp(log_R[b]))`.
Returns a `Makie.Figure`.

Per-draw spaghetti is built as a long-form `DataFrame` and drawn via
AlgebraOfGraphics with `group = :draw`, which is the idiomatic way to
spell "one line per draw" once the data is tidy.

Knot dates default to `BIN_EDGES` (data.jl). Pass `t0` (with optional
`edges`) to recover the calendar dates from a real-time fit whose knot
grid was truncated at the observation cut-off.

# Arguments
- `chn`: a sampled chain from [`joint_model`](@ref).

# Keyword Arguments
- `n_draws_plot`: target number of posterior draws to thin to before
  plotting.
- `ymax`: upper limit of the R(t) y-axis.
- `t0`: time origin (`Date`) for the fit. When supplied, the x-axis is
  mapped to calendar dates `t0 + Day(edges[b])`. Defaults to the package
  retrospective grid `BIN_EDGES`.
- `edges`: knot positions in days from `t0`, as returned by
  [`prepare_rt_edges`](@ref). Defaults to `prepare_rt_edges(t0)` when `t0`
  is supplied.
"""
function plot_rt(chn; n_draws_plot::Int = 100, ymax::Real = 4.0,
        t0::Union{Nothing, Date} = nothing,
        edges::Union{Nothing, AbstractVector} = nothing)
    log_R = vector_chain(chn, :log_R)
    n_draws = length(log_R[1])
    step = max(1, n_draws ÷ n_draws_plot)
    idx = collect(1:step:n_draws)

    knot_dates = _knot_dates(length(log_R); t0, edges)

    # Long form: one row per (draw, knot). `draw` is the only grouping
    # variable, so AoG draws one piecewise-linear trajectory per draw.
    df = DataFrame(
        draw = repeat(idx; inner = length(knot_dates)),
        date = repeat(knot_dates; outer = length(idx)),
        R = vcat([[exp(log_R[b][d]) for b in eachindex(log_R)] for d in idx]...)
    )

    return _with_theme() do
        fig = Figure(; size = (1000, 500))
        spec = data(df) *
               mapping(:date => "Date",
                   :R => "R(t)",
                   group = :draw => nonnumeric) *
               visual(Lines, color = (:steelblue, 0.25), linewidth = 1.6)
        ag = draw!(fig[1, 1], spec;
            axis = (title = "Case reproduction number (weekly knots)",
                limits = (nothing, (0.0, ymax))))
        ax = only(ag).axis
        hlines!(ax, [1.0]; color = :grey, linestyle = :dash)

        _set_date_xticks!(ax, knot_dates)
        fig
    end
end

# Resolve knot calendar dates for a chain with `n_knots` log_R entries.
# Priority: explicit `edges` + `t0` → `prepare_rt_edges(t0)` truncated to
# `n_knots` → `BIN_EDGES[1:n_knots]` (the retrospective default).
function _knot_dates(n_knots::Int; t0::Union{Nothing, Date} = nothing,
        edges::Union{Nothing, AbstractVector} = nothing)
    if edges !== nothing
        t0 === nothing &&
            error("`edges` requires `t0` to convert offsets to dates")
        length(edges) >= n_knots ||
            error("length(edges)=$(length(edges)) < n_knots=$n_knots")
        return Date[t0 + Day(round(Int, edges[b])) for b in 1:n_knots]
    end
    if t0 !== nothing
        e = prepare_rt_edges(t0)
        length(e) >= n_knots ||
            error("default knot grid has $(length(e)) entries; need $n_knots")
        return Date[t0 + Day(round(Int, e[b])) for b in 1:n_knots]
    end
    return BIN_EDGES[1:n_knots]
end

# Pin the x-axis ticks to a readable subset of `knot_dates` rendered as
# ISO date strings. AoG's default Date recipe would print epoch integers
# here because the underlying values are `Date` numerics; switching to an
# explicit tick set keeps the rendered labels human-readable.
function _set_date_xticks!(ax, knot_dates::AbstractVector{Date})
    n = length(knot_dates)
    keep = unique(round.(Int, range(1, n; length = min(n, 7))))
    ax.xticks = (Dates.value.(knot_dates[keep]),
        string.(knot_dates[keep]))
    ax.xticklabelrotation = π / 6
    return ax
end

"""
$(TYPEDSIGNATURES)

Sense-check the per-pair posterior of δ against the fitted population
`Normal(μ_δ, σ_δ)`. For each sourced pair, take the posterior of
`δ_pair = T_inf[secondary] − T_onset[source]` and reduce to its median; then
plot the histogram of those per-pair medians with the population density
overlaid. Returns a `Makie.Figure`.

# Arguments
- `chn`: a sampled chain from [`joint_model`](@ref).
- `d`: the augmented data NamedTuple from [`build_data`](@ref).
"""
function plot_delta_sense_check(chn, d)
    t_inf = vector_chain(chn, :T_inf)
    t_onset = vector_chain(chn, :T_onset)
    μ_δ = _draws(chn, :μ_δ)
    σ_δ = _draws(chn, :σ_δ)

    medians = Float64[]
    for i in 1:d.N
        src = d.source_idx[i]
        src == 0 && continue
        push!(medians, quantile(t_inf[i] .- t_onset[src], 0.5))
    end

    μ_med = quantile(μ_δ, 0.5)
    σ_med = quantile(σ_δ, 0.5)

    df = DataFrame(δ = medians)

    return _with_theme() do
        fig = Figure(; size = (900, 500))
        spec = data(df) *
               mapping(:δ => "δ (days from source onset)") *
               visual(Hist; bins = 15, normalization = :pdf,
                   color = (:steelblue, 0.55),
                   strokecolor = :steelblue, strokewidth = 0.5)
        ag = draw!(fig[1, 1], spec;
            axis = (title = "Per-pair δ vs fitted population Normal",
                ylabel = "density"))
        ax = only(ag).axis
        xs = range(μ_med - 4σ_med, μ_med + 4σ_med; length = 200)
        l_fit = lines!(ax, xs, pdf.(Normal(μ_med, σ_med), xs);
            color = :darkorange, linewidth = 2)
        v_zero = vlines!(ax, [0.0]; color = :grey, linestyle = :dash)
        # Build legend entries manually so the histogram (an AoG layer with
        # no plot handle returned) appears alongside the raw-Makie overlays.
        h_handle = PolyElement(color = (:steelblue, 0.55),
            strokecolor = :steelblue, strokewidth = 0.5)
        Legend(fig[1, 2],
            [h_handle, l_fit, v_zero],
            ["per-pair posterior medians (N = $(length(medians)))",
                "Normal(μ_δ, σ_δ) fitted",
                "source onset"];
            framevisible = false, tellwidth = true)
        fig
    end
end

"""
    plot_z_ppc(model, chn, data; rng = Random.MersenneTwister(1), edges = prepare_rt_edges(data.t0))

Posterior-predictive check for the observed offspring counts `Zobs`. For
each posterior draw `d` and each case `i`, samples a replicated offspring
count via [`posterior_predictive`](@ref) dispatched on
`model.defaults.cases`. For the default [`case_model`](@ref) this is
`Z_rep[i, d] ~ NegativeBinomial(k[d], k[d]/(k[d] + R_i))`, where
`R_i = exp(log_R_at(T_onset[i, d], edges, log_R[:, d]))`.

Joint-draw: `T_onset[i]`, `log_R[:]`, and `k` are taken from the same
posterior draw, so the PPC reflects full posterior uncertainty in case
onset times alongside the time-varying R(t) and dispersion.

Compares the count of cases at each `Z` value (0, 1, 2, …) between the
observed line list and the replicated distribution.

The left panel is a rootogram-style bar chart: bars are observed
frequencies; points + 95% pointwise CrI lines are replicated frequencies.
The right column is three stacked subpanels, one per discrete test
statistic (`sum(Z)`, `max(Z)`, `count(Z == 0)`). Each subpanel shows the
histogram of the replicated statistic and the observed value as a dashed
vertical rule. For numeric values (observed, replicated median + 95%
CrI, two-sided Bayesian posterior-predictive p-value) see the companion
`z_ppc_summary`.

Returns a `Makie.Figure`.
"""
# Joint-draw posterior-predictive replication of Z. Returns an
# `(n_draws × N)` matrix where each row is one replicated line list under
# the model's Negative-Binomial likelihood, using the same draw of
# `(T_inf, log_R, k)`. Dispatches `posterior_predictive` on
# `model.defaults.cases` so alternative obs submodels are honoured.
function _z_ppc_replicate(model, chn, d; rng = Random.MersenneTwister(1),
        edges = nothing)
    k_draws = _draws(chn, :k)
    log_R = vector_chain(chn, :log_R)
    t_onset = vector_chain(chn, :T_onset)
    n_draws = length(k_draws)
    N = d.N
    Zobs = d.Zobs

    knots = edges === nothing ?
            prepare_rt_edges(d.t0)[1:length(log_R)] : edges

    cases_sub = model.defaults.cases

    Z_rep = Matrix{Int}(undef, n_draws, N)
    for d_idx in 1:n_draws
        logR_d = [log_R[b][d_idx] for b in eachindex(log_R)]
        k_d = k_draws[d_idx]
        for i in 1:N
            t_i = t_onset[i][d_idx]
            R_i = exp(log_R_at(t_i, knots, logR_d))
            Z_rep[d_idx, i] = posterior_predictive(cases_sub, rng,
                k_d, Zobs[i], R_i, 1.0, 1.0)
        end
    end
    return Z_rep
end

"""
$(TYPEDSIGNATURES)

Companion to `plot_z_ppc` returning a `DataFrame` of numeric
posterior-predictive summaries for three discrete test statistics —
`sum(Z)`, `max(Z)`, and `count(Z = 0)`. Replicates `Z_rep` jointly in
`(T_inf, log_R, k)` to match `plot_z_ppc`. Columns: `statistic`,
`observed`, `rep_median`, `rep_lower_95`, `rep_upper_95`, `p_ppp`, where
`p_ppp = 2 · min(P(T_rep ≥ T_obs), P(T_rep ≤ T_obs))` is the two-sided
Bayesian posterior-predictive p-value.

# Arguments
- `model`: the [`joint_model`](@ref) instance that produced `chn`; its
  `defaults.cases` field dispatches the per-case posterior-predictive
  draw.
- `chn`: a sampled chain from [`joint_model`](@ref).
- `d`: the augmented data NamedTuple from [`build_data`](@ref).

# Keyword Arguments
- `rng`: RNG used to replicate the offspring counts.
- `edges`: weekly knot edges; defaults to `nothing`, in which case
  the knots are taken from `prepare_rt_edges(d.t0)[1:length(log_R)]`
  inside `_z_ppc_replicate`.
"""
function z_ppc_summary(model, chn, d;
        rng = Random.MersenneTwister(1),
        edges = nothing)
    Z_rep = _z_ppc_replicate(model, chn, d; rng, edges)
    n_draws = size(Z_rep, 1)
    Zobs = d.Zobs

    stats = [
        ("sum(Z)", [sum(view(Z_rep, j, :)) for j in 1:n_draws],
            Float64(sum(Zobs))),
        ("max(Z)", [maximum(view(Z_rep, j, :)) for j in 1:n_draws],
            Float64(maximum(Zobs))),
        ("count(Z = 0)",
            [count(==(0), view(Z_rep, j, :)) for j in 1:n_draws],
            Float64(count(==(0), Zobs)))
    ]

    rows = map(stats) do (name, rep, obs)
        rep_f = Float64.(rep)
        p_ge = mean(rep_f .>= obs)
        p_le = mean(rep_f .<= obs)
        (statistic = name,
            observed = obs,
            rep_median = quantile(rep_f, 0.5),
            rep_lower_95 = quantile(rep_f, 0.025),
            rep_upper_95 = quantile(rep_f, 0.975),
            p_ppp = 2 * min(p_ge, p_le))
    end
    return DataFrame(rows)
end

"""
$(TYPEDSIGNATURES)

Posterior predictive check for offspring counts `Z`. Replicates `Z` from
the joint posterior and overlays the simulated counts against the
observed counts per case.

# Arguments
- `model`: the [`joint_model`](@ref) instance that produced `chn`; its
  `defaults.cases` field dispatches the per-case posterior-predictive
  draw.
- `chn`: a sampled chain from [`joint_model`](@ref).
- `d`: the augmented data NamedTuple from [`build_data`](@ref).

# Keyword Arguments
- `rng`: RNG used for the posterior draws.
- `edges`: weekly knot edges; defaults to `nothing`, in which case
  the knots are taken from `prepare_rt_edges(d.t0)[1:length(log_R)]`
  inside `_z_ppc_replicate`.
"""
function plot_z_ppc(model, chn, d;
        rng = Random.MersenneTwister(1),
        edges = nothing)
    N = d.N
    Z_rep = _z_ppc_replicate(model, chn, d; rng, edges)
    n_draws = size(Z_rep, 1)

    Zobs = d.Zobs
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
    rep_med = [quantile(view(rep_counts, :, j), 0.5)
               for j in eachindex(z_values)]
    rep_lo = [quantile(view(rep_counts, :, j), 0.025)
              for j in eachindex(z_values)]
    rep_hi = [quantile(view(rep_counts, :, j), 0.975)
              for j in eachindex(z_values)]

    # Aggregate test statistics per draw.
    sum_rep = [sum(view(Z_rep, d_idx, :)) for d_idx in 1:n_draws]
    max_rep = [maximum(view(Z_rep, d_idx, :)) for d_idx in 1:n_draws]
    zeros_rep = [count(==(0), view(Z_rep, d_idx, :)) for d_idx in 1:n_draws]
    sum_obs = sum(Zobs)
    max_obs = maximum(Zobs)
    zeros_obs = count(==(0), Zobs)

    return _with_theme() do
        fig = Figure(; size = (1300, 700))

        ax1 = Axis(fig[1, 1];
            xlabel = "Offspring count Z",
            ylabel = "Number of cases",
            title = "Z by value: observed vs replicated",
            xticks = collect(z_values))
        b_obs = barplot!(ax1, collect(z_values), Float64.(obs_counts);
            color = (:steelblue, 0.55),
            strokecolor = :steelblue, strokewidth = 0.5)
        # Error bars for the replicated 95% CrI plus a median marker.
        rangebars!(ax1, collect(z_values), rep_lo, rep_hi;
            color = :darkorange, whiskerwidth = 8)
        s_rep = scatter!(ax1, collect(z_values), rep_med;
            color = :darkorange, markersize = 8)

        # Right column: one stacked subpanel per discrete test statistic.
        # Numeric values (observed, replicated median + 95% CrI, posterior-
        # predictive p) are in `z_ppc_summary` rather than overlaid here.
        right = fig[1, 2] = GridLayout()
        stats = [
            ("sum(Z)", Float64.(sum_rep), Float64(sum_obs), :steelblue),
            ("max(Z)", Float64.(max_rep), Float64(max_obs), :darkorange),
            ("count(Z = 0)", Float64.(zeros_rep),
                Float64(zeros_obs), :seagreen)
        ]
        for (k, (name, rep, obs, colour)) in enumerate(stats)
            ax = Axis(right[k, 1];
                xlabel = k == length(stats) ? "Test statistic value" : "",
                ylabel = "density",
                title = name)
            hist!(ax, rep; bins = 30, normalization = :pdf,
                color = (colour, 0.45))
            vlines!(ax, [obs]; color = colour,
                linestyle = :dash, linewidth = 2)
        end

        Legend(fig[2, 1],
            [b_obs, s_rep],
            ["observed", "replicated (median + 95% CrI)"];
            orientation = :horizontal, framevisible = false,
            tellheight = true, tellwidth = false)
        rowsize!(fig.layout, 2, Auto(0.08))
        # Give the rootogram a bit more horizontal real estate than the
        # narrow stacked stat panels.
        colsize!(fig.layout, 1, Auto(1.4))
        fig
    end
end

"""
$(TYPEDSIGNATURES)

Sense-check the per-case posterior of the incubation period against the
fitted population `LogNormal(μ_inc, σ_inc)`. For each case, takes the
posterior of `inc_i = T_onset[i] − T_inf[i]` and reduces to its median;
plots the histogram of those per-case medians with the median PDF (and
95% pointwise ribbon) of the population LogNormal overlaid. Returns a
`Makie.Figure`.

# Arguments
- `chn`: a sampled chain from [`joint_model`](@ref).
- `data`: the augmented data NamedTuple from [`build_data`](@ref).

# Keyword Arguments
- `n_density_draws`: target number of posterior draws to thin to when
  building the population density ribbon.
"""
function plot_inc_sense_check(chn, data; n_density_draws::Int = 200)
    t_inf = vector_chain(chn, :T_inf)
    t_onset = vector_chain(chn, :T_onset)
    μ_inc = _draws(chn, :μ_inc)
    σ_inc = _draws(chn, :σ_inc)

    medians = [quantile(t_onset[i] .- t_inf[i], 0.5) for i in 1:data.N]

    # Thin posterior draws of (μ_inc, σ_inc) for the density ribbon.
    step = max(1, length(μ_inc) ÷ n_density_draws)
    idx = 1:step:length(μ_inc)
    dists = [LogNormal(μ_inc[i], σ_inc[i]) for i in idx]

    upper = max(maximum(medians) * 1.2, 60.0)
    xs = range(0.5, upper; length = 200)
    med, lo, hi = _density_band(xs, dists)

    return _with_theme() do
        fig = Figure(; size = (900, 500))
        ax = Axis(fig[1, 1];
            xlabel = "Incubation period (days)",
            ylabel = "density",
            title = "Per-case Inc vs fitted population LogNormal")
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
$(TYPEDSIGNATURES)

Faceted overlay of marginal posterior histograms across fits. Designed for
side-by-side comparison of multiple chains: one row per value of `row_col`,
one column per value of `col_col`, with one coloured histogram per `:fit`
level inside each panel.

`df` must be long-form with at least the columns `:value`, `:fit`, plus
the `row_col` and `col_col` columns. The histogram is normalised to a
PDF so overlaid fits with different sample sizes share a y-axis.

Used by the real-time monitoring walkthrough (rows = `:obs_date`,
columns = `:param`) and the prior sensitivity walkthrough (rows =
`:scenario`, columns = `:param`).

# Arguments
- `df`: long-form `DataFrame` of posterior draws.

# Keyword Arguments
- `size_kw`: figure size tuple passed to AoG's `figure` keyword.
- `row_col`: name of the column used for the facet rows. Defaults to
  `:obs_date`.
- `col_col`: name of the column used for the facet columns. Defaults to
  `:param`.
- `bins`: number of histogram bins per panel.
"""
function plot_marginal_overlay(df; size_kw = (1500, 1200),
        row_col::Symbol = :obs_date, col_col::Symbol = :param,
        layout_col::Union{Nothing, Symbol} = nothing,
        bins::Integer = 30)
    base = data(df) *
           visual(Hist; bins = bins, normalization = :pdf, alpha = 0.4)
    spec = layout_col === nothing ?
           base * mapping(:value => "value", color = :fit => "fit",
        row = row_col, col = col_col) :
           base * mapping(:value => "value", color = :fit => "fit",
        layout = layout_col)
    return draw(spec;
        facet = (linkxaxes = :colwise, linkyaxes = :none),
        figure = (; size = size_kw))
end

"""
$(TYPEDSIGNATURES)

Prior-predictive panel: histograms of Inc, δ, and GI/SI drawn from the
package's independent priors on `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`. Returns a
`Makie.Figure`.

Three histograms faceted by quantity is the kind of plot AoG was built
for: one long-form data frame, `mapping(:value, layout = :panel)`,
`visual(Hist)`. Each panel still has its own viewing window so long
tails don't squash the bars; rather than per-facet axis limits, the
input is pre-clipped to the window for each panel.

# Keyword Arguments
- `n`: number of prior draws per quantity.
- `rng`: RNG used to draw the prior samples.
"""
function plot_prior_predictives(; n::Int = 5000,
        rng = Random.MersenneTwister(0))
    μ_inc = rand(rng, Normal(3.0, 0.5), n)
    σ_inc = abs.(rand(rng, Normal(0.0, 0.5), n))
    μ_δ = rand(rng, Normal(0.0, 5.0), n)
    σ_δ = abs.(rand(rng, Normal(0.0, 1.0), n))
    inc_s = [rand(rng, LogNormal(μ_inc[i], σ_inc[i])) for i in 1:n]
    δ_s = [rand(rng, Normal(μ_δ[i], σ_δ[i])) for i in 1:n]
    gi_s = δ_s .+ inc_s

    # One row per (panel, sample). Each panel has its own clip window so
    # the histogram bars aren't dominated by long tails. `panel_idx` keeps
    # facets in the intended order (Inc, δ, GI/SI); `panel` provides the
    # displayed title.
    panels = [
        (idx = 1, name = "Inc (prior)", window = (0.0, 80.0), samples = inc_s),
        (idx = 2, name = "δ (prior)", window = (-25.0, 25.0), samples = δ_s),
        (idx = 3, name = "GI / SI (prior)",
            window = (-30.0, 80.0), samples = gi_s)
    ]
    df = vcat([DataFrame(panel_idx = p.idx,
                   value = filter(x -> p.window[1] <= x <= p.window[2], p.samples))
               for p in panels]...)
    title_pairs = [p.idx => p.name for p in panels]

    return _with_theme() do
        fig = Figure(; size = (1500, 400))
        spec = data(df) *
               mapping(:value => "value",
                   col = :panel_idx => AlgebraOfGraphics.renamer(title_pairs...)) *
               visual(Hist; bins = 100, normalization = :pdf,
                   color = :steelblue)
        draw!(fig[1, 1], spec; facet = (linkxaxes = :none, linkyaxes = :none))
        fig
    end
end

function _save_figure(fig, path)
    Makie.save(path, fig; px_per_unit = 2.0)
    return path
end
