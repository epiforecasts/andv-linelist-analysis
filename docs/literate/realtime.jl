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
# This page checks whether those corrections do their job by overlaying three fits at the same outbreak cut-off date.

# ## The three fits

using Hantavirus
using Dates: Date

obs_date = Date("2018-12-31")
ll       = load_linelist()

ll_truth = filter_by_exposure(ll, obs_date)   # cases known infected by obs_date
ll_naive = filter_realtime(ll, obs_date)      # cases known observed by obs_date

(retained_truth = nrow(ll_truth), retained_naive = nrow(ll_naive))

# | Fit | Linelist | `obs_time` | Interpretation |
# |---|---|---|---|
# | **Counterfactual retro** | `ll_truth` (exposure-filtered) | `nothing` | What an analyst would estimate at obs_date if they had a time machine — full Inc/δ info for every case infected by then. Not realisable in real time. |
# | **Corrected real-time** | `ll_naive` (onset-filtered) | `obs_date` | What `joint_model` offers in real time, with truncation and cluster-completeness corrections active. |
# | **Naive real-time** | `ll_naive` (onset-filtered) | `nothing` | What you'd get if you just filtered the linelist and fit without corrections — the bias the corrections are supposed to fix. |

samples = 300
chains  = 2
seed    = 20260512
tmp     = (output = mktempdir(), figures = mktempdir())

chn_truth, post_truth = analyse(; data = ll_truth,
                                  samples, chains, seed, tmp...,)
chn_rt,    post_rt    = analyse(; data = ll_naive, obs_time = obs_date,
                                  samples, chains, seed, tmp...,)
chn_naive, post_naive = analyse(; data = ll_naive,
                                  samples, chains, seed, tmp...,)
nothing #hide

# ## Population posteriors
#
# For each of `(μ_inc, σ_inc, μ_δ, σ_δ, k)`, the three posteriors are overlaid on the same panel.
# The corrected real-time fit should track the counterfactual retro; the naive fit should diverge from it in the direction of the right-truncation bias (shorter Inc, more negative δ).

using Plots: plot, histogram!, vline!, hline!, plot!, savefig
using Statistics: median

function overlay_marginals(fits, title)
    params = [(:μ_inc, "μ_inc"), (:σ_inc, "σ_inc"),
              (:μ_δ, "μ_δ"),     (:σ_δ, "σ_δ"),
              (:k,   "k")]
    panels = Any[]
    colours = [:steelblue, :darkorange, :crimson]
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
    ("naive real-time",      post_naive),
]
overlay_marginals(fits, "Posterior marginals at $(obs_date)")

# ## R(t) overlay
#
# Posterior medians with 80% CrI ribbons for each bin, the three fits overlaid on a single axis.

function rt_quantiles(post)
    return (lo  = [quantile(exp.(post.log_R_chain[b]), 0.10) for b in eachindex(post.log_R_chain)],
            med = [quantile(exp.(post.log_R_chain[b]), 0.50) for b in eachindex(post.log_R_chain)],
            hi  = [quantile(exp.(post.log_R_chain[b]), 0.90) for b in eachindex(post.log_R_chain)])
end

plt = plot(; xlabel = "Bin index", ylabel = "R(t) (80% CrI)",
             title  = "R(t) at $(obs_date)",
             legend = :topright, ylims = (0.0, 4.0))
colours = [:steelblue, :darkorange, :crimson]
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
# - **Counterfactual retro vs corrected real-time:** these should agree closely if the corrections are doing their job.
#   Any systematic difference is residual bias from approximations (e.g. the cluster-completeness assumes onset-to-report is zero, that source attribution is correct, etc.).
# - **Corrected vs naive real-time:** the gap is the bias the model fixes.
#   The naive fit will tend to underestimate Inc (long incubators not yet observed), pull μ_δ negative (only short / pre-symptomatic δ observed), and underestimate R(t) in the bins closest to the cut-off.

# ## Caveats
#
# The real-time corrections handle three specific biases — long-incubation cases, late transmissions, and incomplete clusters.
# Not corrected:
#
# - geographic / severity / surveillance reporting biases,
# - the onset-to-report delay (only chain completion is modelled),
# - general under-ascertainment,
# - incomplete source attribution.
#
# The counterfactual retrospective fit on this page also uses the eventual exposure attribution to filter cases, which assumes perfect retrospective attribution — in practice some attribution work happens after the cut-off too.
