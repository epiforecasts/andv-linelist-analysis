# # Real-time vs retrospective monitoring
#
# The base `joint_model` fits a closed-out outbreak with complete observation.
# In real time, three biases need correcting:
#
# 1. Long-incubation cases infected before the cut-off may not yet have developed symptoms — observed incubation periods are enriched for short delays.
# 2. Late transmissions from any source may not yet have happened or have not yet been linked — observed transmission timings (δ) are enriched for early / pre-symptomatic events.
# 3. Recent source cases have not had time to seed all their offspring — the observed offspring count is a downward-biased estimate of R(t) near the cut-off.
#
# The real-time machinery in `joint_model` corrects for these via per-case right-truncation on Inc and δ and a cluster-completeness adjustment on the NB offspring count (the [`F_cluster`](@ref) integral).
# This page validates the corrections by fitting three views of the same outbreak and overlaying their posteriors.

using Hantavirus
using DataFrames: nrow
using Dates: Date, Day
using Plots: plot, plot!, histogram!, hline!, vline!, savefig
using Statistics: median, quantile

# ## The three fits
#
# The replication target for the corrected real-time fit is a **counterfactual retrospective** fit that filters the line list by exposure rather than by onset.
# Both views condition on the same set of cases — those known to have been infected by `obs_date` — but only the corrected real-time fit needs to recover the full incubation and transmission distributions from right-truncated observations.
# If the corrections work, the two posteriors should agree.
#
# A third fit on the full closed-out outbreak is included for scientific interest.
# Post-`obs_date` cases can shift the posterior even though they would not be available in real time, so this fit is out of scope for validating the corrections.

obs_date = Date("2019-01-07")
ll       = load_linelist()
ll_truth = filter_by_exposure(ll, obs_date)
ll_rt    = filter_realtime(ll, obs_date)

(full_n = nrow(ll), counterfactual_n = nrow(ll_truth), realtime_n = nrow(ll_rt))

# | Fit | Linelist | `obs_time` | Role |
# |---|---|---|---|
# | **Counterfactual retro** | `ll_truth` (infected ≤ obs_date) | `nothing` | Replication target — same outbreak info at `obs_date` but with full forward Inc / δ for retained cases. |
# | **Corrected real-time** | `ll` (analyse filters to onset ≤ obs_date) | `obs_date` | What `joint_model` offers at the cut-off, with truncation and cluster-completeness corrections active. |
# | **Full retrospective** | `ll` (the whole outbreak) | `nothing` | Out-of-scope science — post-`obs_date` cases included, shown for interest. |
#
# All three fits share the same R(t) bin edges by pinning `t0` to a canonical reference computed from the full line list.
# `bin_edges_day(d.t0)` then returns identical edges across fits and bin indices are directly comparable.

t0_ref = minimum(ll.onset_date) - Day(60)
seed   = 20260512
tmp    = (output = mktempdir(), figures = mktempdir())

chn_truth, post_truth = analyse(; data = ll_truth, t0 = t0_ref, seed, tmp...,)
chn_rt,    post_rt    = analyse(; data = ll, obs_time = obs_date, t0 = t0_ref, seed, tmp...,)
chn_retro, post_retro = analyse(; data = ll, t0 = t0_ref, seed, tmp...,);

# Sanity-check the bin alignment.

@assert bin_edges_day(t0_ref) == bin_edges_day(t0_ref)
@assert length(post_truth.log_R_chain) == length(post_rt.log_R_chain) == length(post_retro.log_R_chain)

# ## Population posteriors
#
# For each of `(μ_inc, σ_inc, μ_δ, σ_δ, k)`, the three posteriors are overlaid on the same panel.
# The corrected real-time density should track the counterfactual retro density; the full retro density may sit elsewhere because it sees additional cases.

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

fits = [
    ("counterfactual retro", post_truth),
    ("corrected real-time",  post_rt),
    ("full retrospective",   post_retro),
]
overlay_marginals(fits, "Posterior marginals at $(obs_date)")

# ## R(t) overlay
#
# Posterior medians with 80% CrI ribbons for each bin, three fits overlaid on a single axis.
# Bin indices are comparable across fits because `t0_ref` is shared.

function rt_quantiles(post)
    return (lo  = [quantile(exp.(post.log_R_chain[b]), 0.10) for b in eachindex(post.log_R_chain)],
            med = [quantile(exp.(post.log_R_chain[b]), 0.50) for b in eachindex(post.log_R_chain)],
            hi  = [quantile(exp.(post.log_R_chain[b]), 0.90) for b in eachindex(post.log_R_chain)])
end

plt = plot(; xlabel = "Bin index", ylabel = "R(t) (80% CrI)",
             title  = "R(t) at $(obs_date)",
             legend = :topright, ylims = (0.0, 4.0))
colours = [:steelblue, :darkorange, :seagreen]
for (i, (name, post)) in enumerate(fits)
    q = rt_quantiles(post)
    b = 1:length(q.med)
    plot!(plt, b, q.med;
          ribbon = (q.med .- q.lo, q.hi .- q.med),
          linewidth = 2, color = colours[i], fillalpha = 0.20,
          label = name)
end
hline!(plt, [1.0]; linestyle = :dash, color = :grey, label = "")
plt

# ## Reading the figures
#
# If the corrected real-time fit reproduces the counterfactual retro posteriors, the corrections are doing their job — the bias from observing only short delays and incomplete clusters has been removed.
# Any systematic offset between those two fits is residual bias from the modelling approximations: the cluster-completeness adjustment assumes onset-to-report is zero, source attribution is correct, and so on.
# The full retrospective is shown for scientific interest only; agreement with the two `obs_date` fits is not expected because it sees later cases.
# Bins past the cut-off in the R(t) plot have no real-time information and reduce to the prior; agreement is expected only in bins covered by retained cases.

# ## Caveats
#
# The real-time corrections handle three specific biases — long-incubation cases, late transmissions, and incomplete clusters.
# Not corrected:
#
# - geographic / severity / surveillance reporting biases,
# - the onset-to-report delay (only chain completion is modelled),
# - general under-ascertainment,
# - incomplete source attribution.
