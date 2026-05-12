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
# This page validates the corrections by fitting the full outbreak retrospectively and the same outbreak at a mid-outbreak cut-off in real-time mode, then overlaying the posteriors.
# If the corrections work, the real-time posteriors at the cut-off should be consistent with the full retrospective posteriors (broader, since the real-time fit sees less data, but centred in the same place).

using Hantavirus
using DataFrames: nrow
using Dates: Date
using Plots: plot, plot!, histogram!, hline!, vline!, savefig
using Statistics: median, quantile

# ## The two fits

obs_date = Date("2018-12-31")
ll       = load_linelist()
ll_rt    = filter_realtime(ll, obs_date)

(retrospective_n = nrow(ll), realtime_n = nrow(ll_rt))

# | Fit | Linelist | `obs_time` | Interpretation |
# |---|---|---|---|
# | **Full retrospective** | `ll` (the whole outbreak) | `nothing` | Gold standard — complete observation, no truncation correction needed. |
# | **Corrected real-time** | `ll_rt` (onset ≤ obs_date) | `obs_date` | What `joint_model` offers at the cut-off, with truncation and cluster-completeness corrections active. |

seed = 20260512
tmp  = (output = mktempdir(), figures = mktempdir())

chn_retro, post_retro = analyse(; data = ll,    seed, tmp...,)
chn_rt,    post_rt    = analyse(; data = ll_rt, obs_time = obs_date, seed, tmp...,);

# ## Population posteriors
#
# For each of `(μ_inc, σ_inc, μ_δ, σ_δ, k)`, the two posteriors are overlaid on the same panel.
# If the real-time corrections are working, the corrected real-time posterior should be centred on the full retrospective posterior; the real-time density will be broader because it conditions on less data.

function overlay_marginals(fits, title)
    params = [(:μ_inc, "μ_inc"), (:σ_inc, "σ_inc"),
              (:μ_δ, "μ_δ"),     (:σ_δ, "σ_δ"),
              (:k,   "k")]
    panels = Any[]
    colours = [:steelblue, :darkorange]
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
    ("full retrospective",  post_retro),
    ("corrected real-time", post_rt),
]
overlay_marginals(fits, "Posterior marginals at $(obs_date)")

# ## R(t) overlay
#
# Posterior medians with 80% CrI ribbons for each bin, the two fits overlaid on a single axis.

function rt_quantiles(post)
    return (lo  = [quantile(exp.(post.log_R_chain[b]), 0.10) for b in eachindex(post.log_R_chain)],
            med = [quantile(exp.(post.log_R_chain[b]), 0.50) for b in eachindex(post.log_R_chain)],
            hi  = [quantile(exp.(post.log_R_chain[b]), 0.90) for b in eachindex(post.log_R_chain)])
end

plt = plot(; xlabel = "Bin index", ylabel = "R(t) (80% CrI)",
             title  = "R(t) at $(obs_date)",
             legend = :topright, ylims = (0.0, 4.0))
colours = [:steelblue, :darkorange]
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
# If the corrected real-time fit produces population posteriors centred on the full retrospective posteriors, the corrections are doing their job — the bias from observing only short delays and incomplete clusters has been removed.
# Any systematic offset is residual bias from the modelling approximations: the cluster-completeness adjustment assumes onset-to-report is zero, source attribution is correct, and so on.
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
