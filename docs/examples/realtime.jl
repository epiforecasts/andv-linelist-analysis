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
# Only once the delay parameters look well-identified does the page move on to the joint fits and their downstream diagnostics, R(t) overlays, population marginals, and controlled-outbreak projection.
#
# Real-time-specific caveats are listed in the [Limitations page](limitations.md)
# under the "Real-time fitting caveats" heading.

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
ll = load_linelist();
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

preps = [_prepare_at(ll, obs_date) for obs_date in obs_dates];

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
end;

# ## Delays-only population posteriors per cut-off
#
# Overlaid marginals of `(μ_inc, σ_inc, μ_δ, σ_δ)` from the delays-only
# fits at each cut-off, one row per `obs_date`. If the corrected
# real-time density tracks the counterfactual retro density on each
# panel, the right-truncation on Inc and δ alone is enough to recover
# the population delays from the truncated observations — and any
# pathology in the corresponding joint fit lives in the R(t) /
# `case_model` half rather than in the delay submodels.

let
    params = [(:μ_inc, "μ_inc"), (:σ_inc, "σ_inc"),
        (:μ_δ, "μ_δ"), (:σ_δ, "σ_δ")]
    colours = [:steelblue, :darkorange]
    fig = Figure(; size = (1500, 700))
    for (r, fit) in enumerate(delays_fits)
        row_fits = [
            ("counterfactual retro", fit.chn_dly_truth),
            ("corrected real-time", fit.chn_dly_rt)
        ]
        for (c, (key, label)) in enumerate(params)
            ax = Axis(fig[r, c];
                xlabel = label, ylabel = "density",
                title = c == 1 ?
                        "obs_date = $(fit.obs_date)" : "")
            for (i, (name, chn)) in enumerate(row_fits)
                draws = vec(collect(chn[key]))
                hist!(ax, draws;
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

# ## Joint fits
#
# Once the delays-only fits are in hand, fit the full joint model:
# one full retrospective fit on the closed-out line list (shared
# comparator across cut-offs) plus, at each cut-off, a counterfactual
# retro and a corrected real-time joint fit.
#
# The full retrospective fit's standalone diagnostics, headline summary
# and data plot are not repeated here — see the analysis walkthrough
# page for those.

d_retro = build_data(ll; t0 = t0_ref);
chn_retro_full = sample_fit(joint_model(d_retro, edges_ref);
    samples = n_samples, chains = n_chains, seed = seed);
post_retro_full = summarise(chn_retro_full);

joint_fits = map(delays_fits) do prep
    @info "Joint at cut-off" prep.obs_date
    chn_truth = sample_fit(joint_model(prep.d_truth, prep.edges_rt);
        samples = n_samples, chains = n_chains, seed = seed)
    chn_rt = sample_fit(joint_model(prep.d_rt, prep.edges_rt);
        samples = n_samples, chains = n_chains, seed = seed)
    merge(prep,
        (; chn_truth, chn_rt,
            post_truth = summarise(chn_truth),
            post_rt = summarise(chn_rt)))
end;

# ## Sampler health across all fits
#
# Combined sampler diagnostics with one row per fit — delays-only and
# joint, retro and real-time, at each cut-off, plus the full retro
# joint reference. R̂ near 1 and zero divergences are the targets.
# If R̂ blows up only on the joint rows at a given cut-off, the delay
# submodels are not to blame and the pathology is in the R(t) /
# `case_model` half of the likelihood.

function _diag_row(chn, obs_date, fit_kind, n_cases)
    d = diagnostics_table(chn)
    return (obs_date = obs_date,
        fit_kind = fit_kind,
        N_cases = n_cases,
        rhat_max = d.rhat_max[1],
        n_divergent = d.divergences[1],
        wall_sec = d.runtime_seconds[1])
end

diag_df = let
    rows = NamedTuple[]
    push!(rows,
        _diag_row(chn_retro_full, missing,
            "full retro (joint)", nrow(ll)))
    for (dly, jnt) in zip(delays_fits, joint_fits)
        push!(rows,
            _diag_row(dly.chn_dly_truth, dly.obs_date,
                "retro (delays only)", dly.n_truth))
        push!(rows,
            _diag_row(dly.chn_dly_rt, dly.obs_date,
                "realtime (delays only)", dly.n_rt))
        push!(rows,
            _diag_row(jnt.chn_truth, jnt.obs_date,
                "retro (joint)", jnt.n_truth))
        push!(rows,
            _diag_row(jnt.chn_rt, jnt.obs_date,
                "realtime (joint)", jnt.n_rt))
    end
    DataFrame(rows)
end

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

# ## Controlled-outbreak projection
#
# At each `obs_date`, two counterfactual predictions for the number of
# future symptomatic cases:
#
# - **Controlled** (`predict_controlled_outbreak`): transmission stops
#   at `obs_time`.
#   Only people already infected (chain started, in incubation) can go
#   on to have onset.
#   Per source `i`, future onsets `~ Poisson(λ_i · (q_i − p_i))` where
#   `q_i = cdf(δ_dist, obs_time − T_onset[i])` is the probability
#   transmission happened by `obs_time` and `p_i =
#   cdf(ConvolvedDelays(inc, δ), obs_time − T_onset[i])` is the
#   probability the chain has completed.
# - **Natural chain** (`predict_natural_chain_outbreak`): current
#   sources keep transmitting at their existing rate, but no second-
#   generation chains form from those new offspring.
#   Per source `i`, future onsets `~ Poisson(λ_i · (1 − p_i))`.
#
# Both share the same Gamma posterior on `λ_i`:
# `λ_i | Z_obs[i], k, R_i, p_i ~ Gamma(k + Z_obs[i], R_i / (k + R_i ·
# p_i))` (scale form).
# The realised count of cases with onset strictly after each cut-off is
# overlaid as a vertical reference; values above the natural-chain band
# imply transmission continued past the cut-off, values below imply it
# stalled.

controlled = map(joint_fits) do fit
    strict = predict_controlled_outbreak(
        fit.chn_rt, fit.post_rt, ll, fit.obs_date, t0_ref)
    natural = predict_natural_chain_outbreak(
        fit.chn_rt, fit.post_rt, ll, fit.obs_date, t0_ref)
    (; fit.obs_date, fit.n_rt,
        strict_samples = strict.future_samples,
        natural_samples = natural.future_samples,
        actual_future = strict.actual_future)
end;

controlled_df = DataFrame(
    obs_date = [c.obs_date for c in controlled],
    n_obs = [c.n_rt for c in controlled],
    actual_future = [c.actual_future for c in controlled],
    strict_median = [Int(round(median(c.strict_samples)))
                     for c in controlled],
    strict_lo10 = [Int(round(quantile(c.strict_samples, 0.10)))
                   for c in controlled],
    strict_hi90 = [Int(round(quantile(c.strict_samples, 0.90)))
                   for c in controlled],
    natural_median = [Int(round(median(c.natural_samples)))
                      for c in controlled],
    natural_lo10 = [Int(round(quantile(c.natural_samples, 0.10)))
                    for c in controlled],
    natural_hi90 = [Int(round(quantile(c.natural_samples, 0.90)))
                    for c in controlled])

#-

let
    fig = Figure(; size = (1500, 400))
    for (j, c) in enumerate(controlled)
        ax = Axis(fig[1, j];
            xlabel = "Future cases", ylabel = "Density",
            title = "obs_date = $(c.obs_date)  (n_obs=$(c.n_rt))")
        hist!(ax, c.strict_samples;
            bins = 30, normalization = :pdf,
            color = (:steelblue, 0.4),
            strokecolor = :steelblue, strokewidth = 1,
            label = "controlled (strict)")
        hist!(ax, c.natural_samples;
            bins = 30, normalization = :pdf,
            color = (:seagreen, 0.3),
            strokecolor = :seagreen, strokewidth = 1,
            label = "natural chain")
        vlines!(ax, [c.actual_future];
            color = :darkorange, linewidth = 3,
            label = "actual = $(c.actual_future)")
        j == 1 && axislegend(ax; position = :rt)
    end
    fig
end

# ### Why the Jan 7 predictive is wider than Dec 31
#
# The Jan 7 cut-off carries more cases but the predictive distribution
# of future onsets is wider, not narrower.
# The mechanism is in the per-source Gamma–Poisson conjugate posterior
# used by `predict_controlled_outbreak`: each observed source `i`
# contributes
# ```
# λ_i | Z_obs[i], k, R_i, p_i  ~  Gamma(k + Z_obs[i],  R_i / (k + R_i · p_i))   (scale form)
# Z_future[i]                  ~  Poisson(λ_i · (1 − p_i))
# ```
# with `p_i = cdf(ConvolvedDelays(inc, δ), obs_time − T_onset[i])` —
# the probability that a chain rooted at source `i` has already
# completed by the cut-off.
# The Gamma scale `R_i / (k + R_i · p_i)` is large for sources with
# small `p_i` (recent onsets near the cut-off): when `p_i ≈ 0`,
# `scale ≈ R_i / k`, and the conditional variance
# `(k + Z_obs[i]) · (R_i / (k + R_i · p_i))²` is correspondingly
# large.
# The Poisson thinning by `(1 − p_i)` then routes essentially all of
# that wide `λ_i` mass into the future window.
# At Jan 7 the line list adds roughly nine secondary cases with onsets
# in the last week before the cut-off (`Δ_i = obs_time − T_onset[i]`
# is just a few days for these): their `p_i` is small, so each
# contributes a high-variance `Gamma · (1 − p_i)` term to the sum,
# inflating the variance of the total despite the larger sample.
# At Dec 31 fewer sources sit in that small-`Δ_i` regime, so the sum
# is dominated by sources with moderate-to-large `p_i` whose scales
# `R_i / (k + R_i · p_i)` are close to `1 / p_i` and whose `(1 − p_i)`
# thinning further damps the contribution to the future window.
# Additionally, the latest-onset sources have few visible offspring,
# so the R(t) posterior near the Jan 7 cut-off is more prior-driven —
# `R_i` is itself wider — which compounds the Gamma scale and the
# resulting Poisson predictive variance.
