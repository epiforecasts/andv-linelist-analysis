# # Real-time vs retrospective monitoring
#
# The base `joint_model` fits a closed-out outbreak with complete observation.
# In real time, three biases need correcting:
#
# 1. Long-incubation cases infected before the cut-off may not yet have developed symptoms — observed incubation periods are enriched for short delays.
# 2. Late transmissions from any source may not yet have happened or have not yet been linked — observed transmission timings (δ) are enriched for early / pre-symptomatic events.
# 3. Recent source cases have not had time to seed all their offspring — the observed offspring count is a downward-biased estimate of R(t) near the cut-off.
#
# The real-time machinery in `joint_model` corrects for these via per-case right-truncation on Inc and δ and an offspring-completeness adjustment on the NB offspring count (the [`F_offspring`](@ref) integral).
# `F_offspring` is the probability that an offspring's `δ + Inc(sec)` chain has completed by the cut-off, conditional on the source's onset time.
# The argument is `obs_time − T_onset[src]`, not `obs_time − T_inf[src]` — the source's own incubation is a sampled latent already scored, so the offspring delay reduces to `δ + Inc(sec)`.
# This page validates the corrections by fitting the same outbreak at three real-time cut-offs and overlaying the resulting R(t) posteriors and population marginals against a counterfactual retrospective and the full closed-out fit.

using TransmissionLinelist
using DataFrames: DataFrame, nrow
using Dates: Date, Day
using Statistics: quantile
using CairoMakie

# ## The three fits per cut-off
#
# At each `obs_date` the replication target for the corrected real-time fit is a **counterfactual retrospective** fit that filters the line list by exposure rather than by onset.
# Both views condition on the same set of cases (those known to have been infected by `obs_date`) but only the corrected real-time fit needs to recover the full incubation and transmission distributions from right-truncated observations.
# If the corrections work, the two posteriors should agree.
#
# A single full retrospective fit on the closed-out outbreak is shared across all three cut-offs and shown for scientific interest.

obs_dates = [Date("2018-12-15"), Date("2018-12-31"), Date("2019-01-07")]
ll        = load_linelist()
t0_ref    = minimum(ll.onset_date) - Day(60)
seed      = 20260512
tmp       = (output = mktempdir(), figures = mktempdir())

# Full retrospective fit (independent of `obs_date`, fit once).

chn_retro, post_retro = analyse(; data = ll, t0 = t0_ref, seed, tmp...,)

# At each cut-off, fit the counterfactual retrospective (exposure-filtered, no `obs_time`) and the corrected real-time view (whole line list, `obs_time = obs_date`).

fits_by_date = map(obs_dates) do obs_date
    ll_truth = filter_by_exposure(ll, obs_date)
    ll_rt    = filter_realtime(ll, obs_date)
    chn_truth, post_truth = analyse(;
        data = ll_truth, t0 = t0_ref, seed, tmp...,)
    chn_rt,    post_rt    = analyse(;
        data = ll, obs_time = obs_date, t0 = t0_ref, seed, tmp...,)
    (; obs_date,
       n_truth = nrow(ll_truth), n_rt = nrow(ll_rt),
       chn_truth, post_truth, chn_rt, post_rt)
end

# All fits share R(t) bin edges by pinning `t0` to `t0_ref` so `bin_edges_day(d.t0)` returns identical edges across fits.

@assert length(post_retro.log_R_chain) ==
        length(fits_by_date[1].post_truth.log_R_chain) ==
        length(fits_by_date[1].post_rt.log_R_chain)

# ## Sampler diagnostics per fit
#
# One row per fit: maximum R̂, divergence count and wall-clock sampling time give a quick read on chain pathology at each cut-off.
# Smaller cut-offs have fewer cases and weaker likelihood, so some divergence and inflated R̂ relative to the full retrospective is expected.

function _diag_row(chn, obs_date, fit_kind, n_cases)
    d = diagnostics_table(chn)
    return (obs_date    = obs_date,
            fit_kind    = fit_kind,
            N_cases     = n_cases,
            rhat_max    = d.rhat_max[1],
            n_divergent = d.divergences[1],
            wall_sec    = d.runtime_seconds[1])
end

diag_rows = NamedTuple[]
push!(diag_rows,
      _diag_row(chn_retro, missing, "full retro", nrow(ll)))
for fit in fits_by_date
    push!(diag_rows,
          _diag_row(fit.chn_truth, fit.obs_date, "retro",
                    fit.n_truth))
    push!(diag_rows,
          _diag_row(fit.chn_rt, fit.obs_date, "realtime",
                    fit.n_rt))
end
diag_df = DataFrame(diag_rows)
show(stdout, MIME"text/plain"(), diag_df)

# ## R(t) per cut-off
#
# One panel per `obs_date` with all three fits overlaid.
# Posterior medians with 80% CrI ribbons; bin indices are comparable across panels because `t0_ref` is shared.

function rt_quantiles(post)
    return (lo  = [quantile(exp.(post.log_R_chain[b]), 0.10)
                   for b in eachindex(post.log_R_chain)],
            med = [quantile(exp.(post.log_R_chain[b]), 0.50)
                   for b in eachindex(post.log_R_chain)],
            hi  = [quantile(exp.(post.log_R_chain[b]), 0.90)
                   for b in eachindex(post.log_R_chain)])
end

let
    colours = [:steelblue, :darkorange, :seagreen]
    fig = Figure(; size = (1500, 500))
    for (j, fit) in enumerate(fits_by_date)
        ax = Axis(fig[1, j];
                  xlabel = "Bin index", ylabel = "R(t) (80% CrI)",
                  title  = "obs_date = $(fit.obs_date)",
                  limits = (nothing, (0.0, 4.0)))
        panel_fits = [
            ("counterfactual retro", fit.post_truth),
            ("corrected real-time",  fit.post_rt),
            ("full retrospective",   post_retro),
        ]
        for (i, (name, post)) in enumerate(panel_fits)
            q = rt_quantiles(post)
            b = collect(1:length(q.med))
            band!(ax, b, q.lo, q.hi; color = (colours[i], 0.2))
            lines!(ax, b, q.med; color = colours[i], linewidth = 2,
                   label = name)
        end
        hlines!(ax, [1.0]; color = :grey, linestyle = :dash)
        j == 1 && axislegend(ax; position = :rt)
    end
    fig
end

# ## Population posteriors per cut-off
#
# Overlaid marginals of `(μ_inc, σ_inc, μ_δ, σ_δ, k)` at each cut-off, one row per `obs_date`.
# If the corrected real-time density tracks the counterfactual retro density on each panel, the corrections are recovering the population from the truncated observations.
# At the earliest cut-offs few cases are retained and the likelihood is weak, so divergence between the two — and from the full retrospective — is expected; agreement should tighten as the cut-off moves later.

let
    params = [(:μ_inc, "μ_inc"), (:σ_inc, "σ_inc"),
              (:μ_δ, "μ_δ"),     (:σ_δ, "σ_δ"),
              (:k,   "k")]
    colours = [:steelblue, :darkorange, :seagreen]

    fig = Figure(; size = (1500, 900))
    for (r, fit) in enumerate(fits_by_date)
        row_fits = [
            ("counterfactual retro", fit.post_truth),
            ("corrected real-time",  fit.post_rt),
            ("full retrospective",   post_retro),
        ]
        for (c, (key, label)) in enumerate(params)
            ax = Axis(fig[r, c];
                      xlabel = label, ylabel = "density",
                      title  = c == 1 ?
                               "obs_date = $(fit.obs_date)" : "")
            for (i, (name, post)) in enumerate(row_fits)
                hist!(ax, getproperty(post, key);
                      bins = 30, normalization = :pdf,
                      color = (colours[i], 0.3),
                      strokecolor = colours[i], strokewidth = 1,
                      label = name)
            end
            r == 1 && c == 1 && axislegend(ax; position = :rt)
        end
    end
    fig
end

# ## Reading the figures
#
# If the corrected real-time fit reproduces the counterfactual retro posteriors at each cut-off, the corrections are doing their job — the bias from observing only short delays and incomplete clusters has been removed.
# Any systematic offset between those two fits is residual bias from the modelling approximations: the offspring-completeness adjustment assumes onset-to-report is zero, source attribution is correct, and so on.
# The full retrospective is shown for scientific interest only; agreement with the cut-off fits is not expected because it sees later cases.
# Bins past the cut-off in the R(t) plot have no real-time information and reduce to the prior; agreement is expected only in bins covered by retained cases.

# ## Caveats
#
# The real-time corrections handle three specific biases — long-incubation cases, late transmissions, and incomplete clusters.
# Not corrected:
#
# - geographic / severity / surveillance reporting biases,
# - the onset-to-report delay (only chain completion is modelled),
# - general under-ascertainment,
# - incomplete source attribution,
# - pre-symptomatic transmission with an unobserved source,
# - ongoing zoonosis.
#
# **Pre-symptomatic transmission with an unobserved source.**
# When `δ < −Inc[src]`, a source's onset can be later than its secondary's onset.
# At an obs_time cut-off the secondary can be in the line list while the source isn't; `filter_realtime` then drops the source attribution and the secondary looks like an apparent index.
# Probably small for ANDV (δ averages near zero with σ_δ ≈ 1) but a real selection effect the current implementation doesn't correct for.
#
# **Ongoing zoonosis.**
# The model treats index (zoonotic) cases as a small starter set; cluster-completeness only thins observed sources, it doesn't add back population members whose Inc hasn't completed yet.
# The current implementation is fine for an outbreak with a few initial zoonotic cases and no ongoing zoonosis (the Epuyén pattern); it would under-count cases if zoonosis continued throughout the outbreak.
