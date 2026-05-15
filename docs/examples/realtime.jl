# # Real-time vs retrospective monitoring
#
# The base `joint_model` fits a closed-out outbreak with complete observation.
# In real time, three biases need correcting:
#
# 1. Long-incubation cases infected before the cut-off may not yet have developed symptoms — observed incubation periods are enriched for short delays.
# 2. Late transmissions from any source may not yet have happened or have not yet been linked — observed transmission timings (δ) are enriched for early / pre-symptomatic events.
# 3. Recent source cases have not had time to seed all their offspring — the observed offspring count is a downward-biased estimate of R(t) near the cut-off.
#
# The real-time machinery in `joint_model` corrects for these via per-case right-truncation on Inc and δ and an offspring-completeness adjustment on the NB offspring count, expressed as the cdf of the [`ConvolvedDelays`](@ref) distribution `δ + Inc(sec)`.
# The adjustment is the probability that an offspring's chain has completed by the cut-off, conditional on the source's onset time.
# The argument is `obs_time − T_onset[src]`, not `obs_time − T_inf[src]` — the source's own incubation is a sampled latent already scored, so the offspring delay reduces to `δ + Inc(sec)`.
# This page validates the corrections by fitting the same outbreak at two real-time cut-offs and overlaying the resulting R(t) posteriors and population marginals against a counterfactual retrospective and the full closed-out fit.
# It also runs a **delays-only diagnostic** at each cut-off — fitting just the incubation and δ submodels — so that if the full joint fit collapses, we can tell whether the pathology lives in the delay submodels or in the R(t) / `case_model` half of the likelihood.
#
# The workflow is sequenced so the delays-only fits at both cut-offs are inspected first.
# If the delay submodels diverge at a cut-off, the joint fit is hopeless and the reader can stop there.
# Only once the delay parameters look well-identified does the page move on to the joint fits and their downstream diagnostics, R(t) overlays, population marginals, posterior-predictive checks, and controlled-outbreak projection.

using TransmissionLinelist
using DataFrames: DataFrame, nrow
using Dates: Dates, Date, Day
using Statistics: quantile, median
using Random
using CairoMakie

Random.seed!(20260512)

# ## Setup
#
# Two cut-offs are used: 31 December 2018 (about three weeks into the
# outbreak) and 7 January 2019 (one week later, by which point the
# joint fit becomes harder to sample).

obs_dates = [Date("2018-12-31"), Date("2019-01-07")]
ll = load_linelist()
t0_ref = minimum(ll.onset_date) - Day(60)
edges_ref = bin_edges_day(t0_ref)
seed = 20260512

# Sampler settings are kept modest (2 chains, 500 samples) to keep the
# documentation build cost low. The package defaults (4 chains, 1000
# samples) are recommended for production use.

n_chains = 2
n_samples = 500

# ## Data preparation per cut-off
#
# For each `obs_date` we build two views of the line list:
#
# - the **counterfactual retrospective** (filtered by exposure: cases
#   known to have been infected by `obs_date`),
# - the **real-time** view (filtered by onset: what an analyst sees at
#   `obs_date`).
#
# Each cut-off also gets its own weekly knot grid `edges_rt` truncated
# at the cut-off so the R(t) posterior does not extend past the
# observation window.

function _prepare_at(ll, obs_date)
    ll_truth = filter_by_exposure(ll, obs_date)
    ll_rt = filter_realtime(ll, obs_date)
    d_truth = build_data(ll_truth; t0 = t0_ref)
    d_rt = build_data(ll_rt; obs_time = obs_date, t0 = t0_ref)
    edges_rt = let
        obs_offset = Float64(Dates.value(obs_date - t0_ref))
        e = edges_ref[edges_ref .<= obs_offset]
        (isempty(e) || e[end] < obs_offset) ? push!(e, obs_offset) : e
    end
    return (; obs_date, ll_truth, ll_rt, d_truth, d_rt, edges_rt,
        n_truth = nrow(ll_truth), n_rt = nrow(ll_rt))
end

preps = [_prepare_at(ll, obs_date) for obs_date in obs_dates]

# ## Delays-only fits
#
# Fit `delays_only_model` on both the counterfactual retrospective and
# the real-time views at each cut-off. These fits drop the R(t) /
# `case_model` NB likelihood and condition on only the incubation and
# δ submodels (plus their truncation in the real-time case). If the
# delay parameters diverge here, the joint fit at the same cut-off has
# no chance.

delays_fits = map(preps) do prep
    @info "Delays-only at cut-off" prep.obs_date
    chn_dly_truth = sample_fit(delays_only_model(prep.d_truth);
        samples = n_samples, chains = n_chains, seed = seed)
    chn_dly_rt = sample_fit(delays_only_model(prep.d_rt);
        samples = n_samples, chains = n_chains, seed = seed)
    merge(prep, (; chn_dly_truth, chn_dly_rt))
end

# ## Delays-only diagnostics
#
# One row per delays-only fit. If R̂ blows up here, the delay submodels
# alone cannot identify the population delays at that cut-off and the
# corresponding joint fit will not save it.

function _diag_row(chn, obs_date, fit_kind, n_cases)
    d = diagnostics_table(chn)
    return (obs_date = obs_date,
        fit_kind = fit_kind,
        N_cases = n_cases,
        rhat_max = d.rhat_max[1],
        n_divergent = d.divergences[1],
        wall_sec = d.runtime_seconds[1])
end

delays_diag_rows = NamedTuple[]
for fit in delays_fits
    push!(delays_diag_rows,
        _diag_row(fit.chn_dly_truth, fit.obs_date,
            "retro (delays only)", fit.n_truth))
    push!(delays_diag_rows,
        _diag_row(fit.chn_dly_rt, fit.obs_date,
            "realtime (delays only)", fit.n_rt))
end
delays_diag_df = DataFrame(delays_diag_rows)

# Posterior medians and 95% CrIs for the four delay parameters across
# the delays-only fits. Stable medians across the retro and real-time
# columns at a given cut-off mean the right-truncation correction is
# recovering the population from the truncated observations.

function _q3(x)
    return (med = quantile(x, 0.5),
        lo = quantile(x, 0.025), hi = quantile(x, 0.975))
end

function _delay_rows(chn, obs_date, fit_kind, n_cases)
    rows = NamedTuple[]
    for name in (:μ_inc, :σ_inc, :μ_δ, :σ_δ)
        q = _q3(vec(collect(chn[name])))
        push!(rows,
            (obs_date = obs_date, fit_kind = fit_kind,
                N_cases = n_cases, param = String(name),
                median = q.med, lo95 = q.lo, hi95 = q.hi))
    end
    return rows
end

delays_param_rows = NamedTuple[]
for fit in delays_fits
    append!(delays_param_rows,
        _delay_rows(fit.chn_dly_truth, fit.obs_date,
            "retro (delays only)", fit.n_truth))
    append!(delays_param_rows,
        _delay_rows(fit.chn_dly_rt, fit.obs_date,
            "realtime (delays only)", fit.n_rt))
end
delays_param_df = DataFrame(delays_param_rows)

# Sense-check panels for the Dec 31 real-time delays-only fit: per-case
# augmented Inc draws and per-pair augmented δ draws against the
# population distributions implied by the posterior. These probe
# whether the delay submodels alone can recover the population delays
# given the augmented latents.

let fit = delays_fits[1]
    plot_inc_sense_check(fit.chn_dly_rt, fit.d_rt)
end

#-

let fit = delays_fits[1]
    plot_delta_sense_check(fit.chn_dly_rt, fit.d_rt)
end

# ## Joint fits
#
# Once the delays-only fits look acceptable, fit the full joint model:
# one full retrospective fit on the closed-out line list (shared
# comparator across cut-offs) plus, at each cut-off, a counterfactual
# retro and a corrected real-time joint fit.
#
# The full retrospective fit's standalone diagnostics, headline summary
# and data plot are not repeated here — see the analysis walkthrough
# page for those.

d_retro = build_data(ll; t0 = t0_ref)
chn_retro_full = sample_fit(joint_model(d_retro, edges_ref);
    samples = n_samples, chains = n_chains, seed = seed)
post_retro_full = summarise(chn_retro_full)

joint_fits = map(delays_fits) do prep
    @info "Joint at cut-off" prep.obs_date
    chn_truth = sample_fit(joint_model(prep.d_truth, edges_ref);
        samples = n_samples, chains = n_chains, seed = seed)
    chn_rt = sample_fit(joint_model(prep.d_rt, prep.edges_rt);
        samples = n_samples, chains = n_chains, seed = seed)
    merge(prep,
        (; chn_truth, chn_rt,
            post_truth = summarise(chn_truth),
            post_rt = summarise(chn_rt)))
end

# ## Joint diagnostics across cut-offs
#
# Combined sampler diagnostics with one row per fit. The delays-only
# rows from the section above are repeated alongside the joint rows so
# that, for each cut-off, the reader can compare R̂ / divergences across
# the two model variants directly. If R̂ blows up only on the joint
# rows, the delay submodels are not to blame and the pathology is in
# the R(t) / `case_model` half of the likelihood.

diag_rows = NamedTuple[]
push!(diag_rows,
    _diag_row(chn_retro_full, missing, "full retro (joint)", nrow(ll)))
for (dly, jnt) in zip(delays_fits, joint_fits)
    push!(diag_rows,
        _diag_row(dly.chn_dly_truth, dly.obs_date,
            "retro (delays only)", dly.n_truth))
    push!(diag_rows,
        _diag_row(dly.chn_dly_rt, dly.obs_date,
            "realtime (delays only)", dly.n_rt))
    push!(diag_rows,
        _diag_row(jnt.chn_truth, jnt.obs_date,
            "retro (joint)", jnt.n_truth))
    push!(diag_rows,
        _diag_row(jnt.chn_rt, jnt.obs_date,
            "realtime (joint)", jnt.n_rt))
end
diag_df = DataFrame(diag_rows)

# Cross-table comparing delay-parameter posterior medians for the
# delays-only and joint fits at each cut-off. Stable medians across
# the two model variants mean adding the R(t) / `case_model` likelihood
# does not move the delay parameters; large moves point to
# identifiability problems in the joint fit.

delay_rows = NamedTuple[]
for (dly, jnt) in zip(delays_fits, joint_fits)
    append!(delay_rows,
        _delay_rows(dly.chn_dly_truth, dly.obs_date,
            "retro (delays only)", dly.n_truth))
    append!(delay_rows,
        _delay_rows(dly.chn_dly_rt, dly.obs_date,
            "realtime (delays only)", dly.n_rt))
    append!(delay_rows,
        _delay_rows(jnt.chn_truth, jnt.obs_date,
            "retro (joint)", jnt.n_truth))
    append!(delay_rows,
        _delay_rows(jnt.chn_rt, jnt.obs_date,
            "realtime (joint)", jnt.n_rt))
end
delay_df = DataFrame(delay_rows)

# ## R(t) per cut-off
#
# One panel per `obs_date` with the corrected real-time fit, the
# counterfactual retrospective, and the full retrospective overlaid.
# Posterior medians with 80% CrI ribbons. Bin indices are comparable
# across panels because `t0_ref` is shared.

function rt_quantiles(post)
    return (
        lo = [quantile(exp.(post.log_R_chain[b]), 0.10)
              for b in eachindex(post.log_R_chain)],
        med = [quantile(exp.(post.log_R_chain[b]), 0.50)
               for b in eachindex(post.log_R_chain)],
        hi = [quantile(exp.(post.log_R_chain[b]), 0.90)
              for b in eachindex(post.log_R_chain)])
end

let
    colours = [:steelblue, :darkorange, :seagreen]
    fig = Figure(; size = (1500, 500))
    for (j, fit) in enumerate(joint_fits)
        ax = Axis(fig[1, j];
            xlabel = "Bin index", ylabel = "R(t) (80% CrI)",
            title = "obs_date = $(fit.obs_date)",
            limits = (nothing, (0.0, 4.0)))
        panel_fits = [
            ("counterfactual retro", fit.post_truth),
            ("corrected real-time", fit.post_rt),
            ("full retrospective", post_retro_full)
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
# Overlaid marginals of `(μ_inc, σ_inc, μ_δ, σ_δ, k)` at each cut-off,
# one row per `obs_date`. If the corrected real-time density tracks
# the counterfactual retro density on each panel, the corrections are
# recovering the population from the truncated observations.

let
    params = [(:μ_inc, "μ_inc"), (:σ_inc, "σ_inc"),
        (:μ_δ, "μ_δ"), (:σ_δ, "σ_δ"),
        (:k, "k")]
    colours = [:steelblue, :darkorange, :seagreen]
    fig = Figure(; size = (1500, 900))
    for (r, fit) in enumerate(joint_fits)
        row_fits = [
            ("counterfactual retro", fit.post_truth),
            ("corrected real-time", fit.post_rt),
            ("full retrospective", post_retro_full)
        ]
        for (c, (key, label)) in enumerate(params)
            ax = Axis(fig[r, c];
                xlabel = label, ylabel = "density",
                title = c == 1 ?
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

# ## Z PPC and pair plot
#
# Illustrative joint-fit diagnostics for the Dec 31 cut-off. The Inc
# and δ sense-check panels under the joint fit (with the R(t) /
# `case_model` likelihood layered on) for comparison against the
# delays-only versions above:

let fit = joint_fits[1]
    plot_inc_sense_check(fit.chn_rt, fit.d_rt)
end

#-

let fit = joint_fits[1]
    plot_delta_sense_check(fit.chn_rt, fit.d_rt)
end

# Pair plot for the population scalars in the Dec 31 joint real-time
# fit. Strong banana-shaped correlations between `(μ_inc, σ_inc)` and
# `k` are typical of the joint likelihood under truncation.

let fit = joint_fits[1]
    plot_pair(fit.chn_rt)
end

# Offspring posterior-predictive check for the Dec 31 joint real-time
# fit (joint fit only — the delays-only fit has no NB likelihood).

let fit = joint_fits[1]
    plot_z_ppc(fit.chn_rt, fit.d_rt)
end

#-

let fit = joint_fits[1]
    z_ppc_summary(fit.chn_rt, fit.d_rt)
end

# ## Controlled-outbreak projection
#
# At each `obs_date`, conditional on the corrected real-time posterior,
# the number of future symptomatic cases assuming **no further
# transmission** after the cut-off.
# This is the "what if control is achieved now" forecast: secondaries
# already infected by `obs_time` continue to complete their incubation
# and become symptomatic, but no new infections occur.
# Each observed source contributes `Z_future[i] ~ Poisson(λ_i (1 − F_off_i))`
# where `λ_i | Z_obs[i]` follows the conjugate Gamma posterior of the
# NB-binomial-thinning model, sharpening the prediction by conditioning
# on the source's already-observed offspring count.
# The realised count of cases with onset strictly after each cut-off is
# overlaid as a vertical reference; values above the predicted band
# imply transmission continued past the cut-off in the actual outbreak,
# values below imply it stalled.

controlled = map(joint_fits) do fit
    res = predict_controlled_outbreak(
        fit.chn_rt, fit.post_rt, ll, fit.obs_date, t0_ref)
    (; fit.obs_date, fit.n_rt, res.future_samples, res.actual_future)
end

controlled_df = DataFrame(
    obs_date = [c.obs_date for c in controlled],
    n_obs = [c.n_rt for c in controlled],
    actual_future = [c.actual_future for c in controlled],
    pred_median = [Int(round(median(c.future_samples)))
                   for c in controlled],
    pred_lo10 = [Int(round(quantile(c.future_samples, 0.10)))
                 for c in controlled],
    pred_hi90 = [Int(round(quantile(c.future_samples, 0.90)))
                 for c in controlled])

#-

let
    fig = Figure(; size = (1500, 400))
    for (j, c) in enumerate(controlled)
        ax = Axis(fig[1, j];
            xlabel = "Future cases (controlled counterfactual)",
            ylabel = "Density",
            title = "obs_date = $(c.obs_date)  (n_obs=$(c.n_rt))")
        hist!(ax, c.future_samples;
            bins = 30, normalization = :pdf,
            color = (:steelblue, 0.4),
            strokecolor = :steelblue, strokewidth = 1)
        vlines!(ax, [c.actual_future];
            color = :darkorange, linewidth = 3,
            label = "actual = $(c.actual_future)")
        axislegend(ax; position = :rt)
    end
    fig
end

# ## Reading the figures
#
# If the corrected real-time fit reproduces the counterfactual retro
# posteriors at each cut-off, the corrections are doing their job —
# the bias from observing only short delays and incomplete clusters
# has been removed.
# Any systematic offset between those two fits is residual bias from
# the modelling approximations: the offspring-completeness adjustment
# assumes onset-to-report is zero, source attribution is correct, and
# so on.
# The full retrospective is shown for scientific interest only;
# agreement with the cut-off fits is not expected because it sees
# later cases.
# Bins past the cut-off in the R(t) plot have no real-time information
# and reduce to the prior; agreement is expected only in bins covered
# by retained cases.
# The delays-only diagnostic is the place to look first if the joint
# fit looks broken at a given cut-off: if the delay parameters agree
# between the delays-only and joint columns at that cut-off, the
# R(t) / `case_model` half is where the pathology lives.

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
