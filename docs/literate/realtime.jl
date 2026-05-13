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
# The source's own incubation period is *not* marginalised here: the source is observed, so `T_onset[src]` is pinned by the sampled latents in the same model evaluation and the offspring delay reduces to `δ + Inc(sec)`.
# The argument is therefore `obs_time − T_onset[src]`, not `obs_time − T_inf[src]`.
# This page validates the corrections by fitting the same outbreak at three real-time cut-offs and overlaying the resulting posteriors against a counterfactual retrospective and the full closed-out fit.

using Hantavirus
using DataFrames: nrow
using Dates: Date, Day
using Plots: plot, plot!, histogram!, hline!, vline!, savefig
using Statistics: median, quantile

# ## The three fits per cut-off
#
# At each `obs_date` the replication target for the corrected real-time fit is a **counterfactual retrospective** fit that filters the line list by exposure rather than by onset.
# Both views condition on the same set of cases — those known to have been infected by `obs_date` — but only the corrected real-time fit needs to recover the full incubation and transmission distributions from right-truncated observations.
# If the corrections work, the two posteriors should agree.
#
# A single full retrospective fit on the closed-out outbreak is shared across all three cut-offs and shown for scientific interest.
# Post-`obs_date` cases can shift the posterior even though they would not be available in real time, so this fit is out of scope for validating the corrections.

obs_dates = [Date("2018-12-15"), Date("2018-12-31"), Date("2019-01-07")]
ll        = load_linelist()
t0_ref    = minimum(ll.onset_date) - Day(60)
seed      = 20260512
tmp       = (output = mktempdir(), figures = mktempdir())

# The full retrospective fit is independent of `obs_date` and is fit once.

chn_retro, post_retro = analyse(; data = ll, t0 = t0_ref, seed, tmp...,)

# For each cut-off, fit the counterfactual retrospective (filter by exposure, no `obs_time`) and the corrected real-time view (whole line list, `obs_time = obs_date`; `analyse` filters internally).

fits_by_date = map(obs_dates) do obs_date
    ll_truth  = filter_by_exposure(ll, obs_date)
    chn_truth, post_truth = analyse(; data = ll_truth, t0 = t0_ref, seed, tmp...,)
    chn_rt,    post_rt    = analyse(; data = ll, obs_time = obs_date, t0 = t0_ref, seed, tmp...,)
    (; obs_date, post_truth, post_rt)
end

# | Fit | Linelist | `obs_time` | Role |
# |---|---|---|---|
# | **Counterfactual retro** | `filter_by_exposure(ll, obs_date)` | `nothing` | Replication target — same outbreak info at `obs_date` but with full forward Inc / δ for retained cases. |
# | **Corrected real-time** | `ll` (analyse filters to onset ≤ obs_date) | `obs_date` | What `joint_model` offers at the cut-off, with truncation and offspring-completeness corrections active. |
# | **Full retrospective** | `ll` (the whole outbreak) | `nothing` | Out-of-scope science — post-`obs_date` cases included, shown for interest; independent of `obs_date` and shared across panels. |
#
# All fits share the same R(t) bin edges by pinning `t0` to `t0_ref` computed from the full line list.
# `bin_edges_day(d.t0)` then returns identical edges across fits and bin indices are directly comparable.

@assert length(post_retro.log_R_chain) ==
        length(fits_by_date[1].post_truth.log_R_chain) ==
        length(fits_by_date[1].post_rt.log_R_chain)

# ## Population posteriors
#
# Population marginals are shown for the latest cut-off only (`2019-01-07`).
# The three earlier cut-offs differ mainly in how much of the outbreak is visible, which is best read from the R(t) facets below.

function overlay_marginals(fits, title)
    params = [(:μ_inc, "μ_inc"), (:σ_inc, "σ_inc"),
              (:μ_δ, "μ_δ"),     (:σ_δ, "σ_δ"),
              (:k,   "k")]
    panels = Any[]
    colours = [:steelblue, :darkorange, :seagreen]
    for (key, label) in params
        p = plot(; xlabel = label, ylabel = "density",
                   legend = (key == :μ_inc ? :topright : false))
        for (i, (name, post)) in enumerate(fits)
            histogram!(p, getproperty(post, key);
                       bins = 30, normalize = :pdf,
                       linecolor = colours[i], fillcolor = colours[i],
                       fillalpha = 0.30, linealpha = 0.9, label = name)
        end
        push!(panels, p)
    end
    return plot(panels...; layout = (2, 3), size = (1200, 700),
                plot_title = title)
end

latest = fits_by_date[end]
fits_latest = [
    ("counterfactual retro", latest.post_truth),
    ("corrected real-time",  latest.post_rt),
    ("full retrospective",   post_retro),
]
overlay_marginals(fits_latest, "Posterior marginals at $(latest.obs_date)")

# ## R(t) per cut-off
#
# A 3-panel layout, one per `obs_date`, with all three fits overlaid in each panel.
# Posterior medians with 80% CrI ribbons for each bin; bin indices are comparable across fits and panels because `t0_ref` is shared.

function rt_quantiles(post)
    return (lo  = [quantile(exp.(post.log_R_chain[b]), 0.10) for b in eachindex(post.log_R_chain)],
            med = [quantile(exp.(post.log_R_chain[b]), 0.50) for b in eachindex(post.log_R_chain)],
            hi  = [quantile(exp.(post.log_R_chain[b]), 0.90) for b in eachindex(post.log_R_chain)])
end

colours = [:steelblue, :darkorange, :seagreen]
panels = Any[]
for (j, fit) in enumerate(fits_by_date)
    panel_fits = [
        ("counterfactual retro", fit.post_truth),
        ("corrected real-time",  fit.post_rt),
        ("full retrospective",   post_retro),
    ]
    p = plot(; xlabel = "Bin index", ylabel = "R(t) (80% CrI)",
               title  = "obs_date = $(fit.obs_date)",
               legend = (j == 1 ? :topright : false),
               ylims  = (0.0, 4.0))
    for (i, (name, post)) in enumerate(panel_fits)
        q = rt_quantiles(post)
        b = 1:length(q.med)
        plot!(p, b, q.med;
              ribbon = (q.med .- q.lo, q.hi .- q.med),
              linewidth = 2, color = colours[i], fillalpha = 0.20,
              label = name)
    end
    hline!(p, [1.0]; linestyle = :dash, color = :grey, label = "")
    push!(panels, p)
end
plot(panels...; layout = (1, 3), size = (1500, 500),
     plot_title = "R(t) across real-time cut-offs")

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
# This is probably small for ANDV (δ averages near zero with σ_δ ≈ 1) but is a real selection effect that the current implementation does not correct for.
#
# **Ongoing zoonosis.**
# The model treats index (zoonotic) cases as a small starter set; cluster-completeness only thins observed sources, it doesn't add back population members whose Inc hasn't completed yet.
# The current implementation is fine for an outbreak with a few initial zoonotic cases and no ongoing zoonosis (the Epuyén pattern); it would under-count cases if zoonosis continued throughout the outbreak.
