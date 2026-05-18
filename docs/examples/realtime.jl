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
using AlgebraOfGraphics: data, mapping, visual, draw
using DataFrames: DataFrame, nrow
using Dates: Dates, Date, Day
using Random
using Statistics: quantile
using CairoMakie
using CairoMakie: Band, Lines, Hist, VLines
using Logging: Logging

Random.seed!(20260512)
# Silence the NUTS "Found initial step size" Info logs that would
# otherwise clutter the rendered example output.
Logging.disable_logging(Logging.Info)

# ## Setup
#
# Two cut-offs are used: 31 December 2018 (about three weeks into the
# outbreak) and 7 January 2019.

obs_dates = [Date("2018-12-31"), Date("2019-01-07")]
ll = load_linelist();
t0_ref = minimum(ll.onset_date) - Day(60)
edges_ref = prepare_rt_edges(t0_ref)
seed = 20260512

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

function prepare_at(ll, obs_date)
    ll_truth = filter_by_exposure(ll, obs_date)
    ll_rt = filter_realtime(ll, obs_date)
    d_truth = build_data(ll_truth; t0 = t0_ref)
    d_rt = build_data(ll_rt; obs_time = obs_date, t0 = t0_ref)
    edges_rt = prepare_rt_edges(t0_ref; obs_time = obs_date)
    return (; obs_date, ll_truth, ll_rt, d_truth, d_rt, edges_rt,
        n_truth = nrow(ll_truth), n_rt = nrow(ll_rt))
end

preps = [prepare_at(ll, obs_date) for obs_date in obs_dates];

# ## Delays-only fits
#
# Fit `delays_only_model` on both the counterfactual retrospective and
# the real-time views at each cut-off. These fits drop the R(t) /
# `case_model` NB likelihood and condition on only the incubation and
# δ submodels (plus their truncation in the real-time case). If the
# delay parameters diverge here, the joint fit at the same cut-off has
# no chance.

delays_fits = map(preps) do prep
    chn_dly_truth = sample_fit(delays_only_model(prep.d_truth);
        seed = seed)
    chn_dly_rt = sample_fit(delays_only_model(prep.d_rt);
        seed = seed)
    merge(prep, (; chn_dly_truth, chn_dly_rt))
end;

# ## Delays-only population posteriors per cut-off
#
# Overlaid marginals of `(μ_inc, σ_inc, μ_δ, σ_δ)` from the delays-only
# fits at each cut-off, one row per `obs_date`. If the corrected
# real-time density tracks the counterfactual retro density on each
# panel, the right-truncation on Inc and δ alone is enough to recover
# the population delays from the truncated observations.

# Shared marginal-overlay helper used by both this delays-only plot
# and the joint-fit version further down. AlgebraOfGraphics spec with
# per-column free x-axis so each parameter's scale is its own.
function plot_marginal_overlay(df; size_kw = (1500, 700))
    spec = data(df) *
           mapping(:value => "value", color = :fit => "fit",
               row = :obs_date, col = :param) *
           visual(Hist; bins = 30, normalization = :pdf, alpha = 0.4)
    return draw(spec;
        facet = (linkxaxes = :colwise, linkyaxes = :none),
        figure = (; size = size_kw))
end

let
    params = [:μ_inc, :σ_inc, :μ_δ, :σ_δ]
    rows = NamedTuple[]
    for fit in delays_fits
        row_fits = [
            ("counterfactual retro", fit.chn_dly_truth),
            ("corrected real-time", fit.chn_dly_rt)
        ]
        for (name, chn) in row_fits, p in params

            for v in vec(collect(chn[p]))
                push!(rows,
                    (obs_date = string(fit.obs_date),
                        fit = name, param = String(p), value = v))
            end
        end
    end
    plot_marginal_overlay(DataFrame(rows))
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
    seed = seed);
post_retro_full = summarise(chn_retro_full);

joint_fits = map(delays_fits) do prep
    m_truth = joint_model(prep.d_truth, prep.edges_rt)
    m_rt = joint_model(prep.d_rt, prep.edges_rt)
    chn_truth = sample_fit(m_truth; seed = seed)
    chn_rt = sample_fit(m_rt; seed = seed)
    merge(prep,
        (; m_truth, m_rt, chn_truth, chn_rt,
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

function diag_row(chn, obs_date, fit_kind, n_cases)
    d = first(diagnostics_table(chn))
    return merge((; obs_date, fit_kind, N_cases = n_cases),
        (; rhat_max = d.rhat_max,
            n_divergent = d.divergences,
            wall_sec = d.runtime_seconds))
end

diag_df = let
    rows = NamedTuple[]
    push!(rows,
        diag_row(chn_retro_full, missing,
            "full retro (joint)", nrow(ll)))
    for (dly, jnt) in zip(delays_fits, joint_fits)
        push!(rows,
            diag_row(dly.chn_dly_truth, dly.obs_date,
                "retro (delays only)", dly.n_truth))
        push!(rows,
            diag_row(dly.chn_dly_rt, dly.obs_date,
                "realtime (delays only)", dly.n_rt))
        push!(rows,
            diag_row(jnt.chn_truth, jnt.obs_date,
                "retro (joint)", jnt.n_truth))
        push!(rows,
            diag_row(jnt.chn_rt, jnt.obs_date,
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

let
    band_rows = DataFrame[]
    for fit in joint_fits
        panel_fits = [
            ("counterfactual retro", fit.post_truth),
            ("corrected real-time", fit.post_rt),
            ("full retrospective", post_retro_full)
        ]
        for (name, post) in panel_fits
            tbl = rt_band(post)
            tbl.fit .= name
            tbl.obs_date .= string(fit.obs_date)
            push!(band_rows, tbl)
        end
    end
    df = reduce(vcat, band_rows)
    band_spec = data(df) *
                mapping(:bin => "Bin index",
                    :lo => "R(t) (80% CrI)", :hi;
                    color = :fit, col = :obs_date) *
                visual(Band; alpha = 0.2)
    line_spec = data(df) *
                mapping(:bin => "Bin index",
                    :med => "R(t) (80% CrI)";
                    color = :fit, col = :obs_date) *
                visual(Lines; linewidth = 2)
    draw(band_spec + line_spec;
        axis = (; limits = (nothing, (0.0, 4.0))),
        figure = (; size = (1500, 500)))
end

# ## Population posteriors per cut-off
#
# Overlaid marginals of `(μ_inc, σ_inc, μ_δ, σ_δ, k)` at each cut-off,
# one row per `obs_date`. If the corrected real-time density tracks
# the counterfactual retro density on each panel, the corrections are
# recovering the population from the truncated observations.

let
    params = [:μ_inc, :σ_inc, :μ_δ, :σ_δ, :k]
    ## Cap `k` at its 99% quantile (pooled across fits) so the long
    ## right tail doesn't compress the other panels visually.
    all_k = Float64[]
    for fit in joint_fits,
        post in
        (fit.post_truth, fit.post_rt, post_retro_full)

        append!(all_k, getproperty(post, :k))
    end
    k_cap = quantile(filter(!isnan, all_k), 0.99)
    rows = NamedTuple[]
    for fit in joint_fits
        row_fits = [
            ("counterfactual retro", fit.post_truth),
            ("corrected real-time", fit.post_rt),
            ("full retrospective", post_retro_full)
        ]
        for (name, post) in row_fits, p in params

            for v in getproperty(post, p)
                p === :k && v > k_cap && continue
                push!(rows,
                    (obs_date = string(fit.obs_date),
                        fit = name, param = String(p), value = v))
            end
        end
    end
    plot_marginal_overlay(DataFrame(rows); size_kw = (1500, 900))
end

# ## Controlled-outbreak projection
#
# At each `obs_date`, three counterfactual predictions for the number
# of future symptomatic cases, each conditional on an *assumed*
# intervention rule. None of these scenarios is a factual claim about
# what happened in Epuyén; they describe what the fitted model implies
# *if* transmission had been interrupted in one of three ways.
#
# - **Strict (assumed intervention at 2018-12-31)**
#   ([`predict_controlled_outbreak`](@ref) with
#   `intervention_time = cutoff_18`): the counterfactual assumes all
#   transmission had stopped on `cutoff_18 = 2018-12-31` (the onset
#   date of the 18th case). Post-cutoff sources contribute zero by
#   construction (`Δ_q[i] < 0` ⇒ `q_i ≈ 0`).
# - **Controlled at obs** ([`predict_controlled_outbreak`](@ref) with
#   the default `intervention_time = nothing`): assumes transmission
#   stops at each cut-off, so `Δ_q = Δ_p` and `π_i = q_i − p_i`.
# - **Natural chain** ([`predict_natural_chain_outbreak`](@ref)):
#   current sources keep transmitting at their existing rate but no
#   second-generation chains form from those new offspring.
#
# See the [Model page](model.md#Real-time-predictions) for the
# Gamma–Poisson conjugate posterior these share and the per-source
# thinning probabilities that distinguish them.
# The realised count of cases with onset strictly after each cut-off
# is overlaid as a vertical reference. It is one realisation drawn
# from the true (unknown) process; the predictive distribution is the
# posterior over what could have happened consistent with the fitted
# model and the assumed intervention rule, so we expect spread around
# the realised count rather than agreement on a point.
#
# ### Per-source intervention rules
#
# The scenarios above set a single intervention date (or none) that
# applies to every source. The `intervention_time` keyword also
# accepts a `Vector{Date}` of length `N`, one cut-off per observed
# source. This is useful when isolation policies differ across cases.
# If we believed cases were being isolated on or shortly after their
# own symptom onset, for example under an aggressive contact-tracing
# policy, we could encode that scenario by passing a per-source
# `Vector{Date}` derived from `ll_rt.onset_date`. We do not have
# direct evidence of per-case isolation in the Epuyén outbreak so we
# do not run this scenario here, but the API supports it.

# `predict_*_outbreak` is pure in the fit (model + chain + posterior +
# the `d` it was fit on). The realised future count comes from a
# separate call to `realised_future_count(ll, obs_date)`, so the
# comparator is decoupled from the prediction.

cutoff_18 = Date("2018-12-31")

controlled = map(joint_fits) do fit
    strict = predict_controlled_outbreak(
        fit.m_rt, fit.chn_rt, fit.post_rt, fit.d_rt;
        obs_time = fit.obs_date, t0 = t0_ref,
        intervention_time = cutoff_18)
    at_obs = predict_controlled_outbreak(
        fit.m_rt, fit.chn_rt, fit.post_rt, fit.d_rt;
        obs_time = fit.obs_date, t0 = t0_ref)
    natural = predict_natural_chain_outbreak(
        fit.m_rt, fit.chn_rt, fit.post_rt, fit.d_rt;
        obs_time = fit.obs_date, t0 = t0_ref)
    (; fit.obs_date, fit.n_rt,
        strict_samples = strict.future_samples,
        at_obs_samples = at_obs.future_samples,
        natural_samples = natural.future_samples,
        at_obs_per_source = at_obs,
        actual_future = realised_future_count(ll, fit.obs_date))
end;

# Per-source breakdown for the controlled-at-obs counterfactual at the
# latest cut-off: shows which sources drive the projected future onset
# count. One row per observed source case in the order they appear in
# the predictor's `d` argument.

per_source_predictive_summary(controlled[end].at_obs_per_source)

controlled_df = let
    rows = map(controlled) do c
        s = summarise_predictive(c.strict_samples)
        a = summarise_predictive(c.at_obs_samples)
        n = summarise_predictive(c.natural_samples)
        (obs_date = c.obs_date,
            n_obs = c.n_rt,
            actual_future = c.actual_future,
            strict_median = s.med,
            strict_lo10 = s.lo,
            strict_hi90 = s.hi,
            at_obs_median = a.med,
            at_obs_lo10 = a.lo,
            at_obs_hi90 = a.hi,
            natural_median = n.med,
            natural_lo10 = n.lo,
            natural_hi90 = n.hi)
    end
    DataFrame(rows)
end

#-

let
    hist_rows = NamedTuple[]
    vline_rows = NamedTuple[]
    for c in controlled
        panel = "obs_date = $(c.obs_date)  (n_obs=$(c.n_rt))"
        ## Cap per-panel x-axis at the 99% quantile of all three
        ## variants pooled, so the long natural-chain tail doesn't
        ## compress the strict and at-obs bulks.
        pooled = vcat(Float64.(c.strict_samples),
            Float64.(c.at_obs_samples),
            Float64.(c.natural_samples))
        cap = quantile(pooled, 0.99)
        for v in c.strict_samples
            v > cap && continue
            push!(hist_rows,
                (panel = panel,
                    kind = "strict (intervention at $(cutoff_18))",
                    value = v))
        end
        for v in c.at_obs_samples
            v > cap && continue
            push!(hist_rows,
                (panel = panel, kind = "controlled at obs", value = v))
        end
        for v in c.natural_samples
            v > cap && continue
            push!(hist_rows,
                (panel = panel, kind = "natural chain", value = v))
        end
        push!(vline_rows,
            (panel = panel,
                kind = "actual = $(c.actual_future)",
                value = Float64(c.actual_future)))
    end
    df_hist = DataFrame(hist_rows)
    df_vline = DataFrame(vline_rows)
    hist_spec = data(df_hist) *
                mapping(:value => "Future cases";
                    color = :kind => "kind", col = :panel) *
                visual(Hist; bins = 30, normalization = :pdf, alpha = 0.4)
    vline_spec = data(df_vline) *
                 mapping(:value; color = :kind => "kind",
                     col = :panel) *
                 visual(VLines; linewidth = 3)
    draw(hist_spec + vline_spec;
        facet = (linkxaxes = :none, linkyaxes = :none),
        figure = (; size = (1500, 400)))
end

# ### Why the Jan 7 predictive is wider than Dec 31
#
# More data does not always mean a tighter future-cases prediction.
# At Jan 7 the line list adds secondary cases with onsets in the last
# week before the cut-off, so for those sources the chain is far from
# complete and most of their expected offspring still lies in the
# future.
# Their R(t) is also more prior-driven (few visible offspring near the
# cut-off), so the per-source rate is itself wide.
# Both effects funnel a lot of variance into the predicted total.
# At Dec 31 the late-onset sources are fewer and most of their
# expected offspring chains have already completed, so each source
# contributes far less to the future window.
#
# ### Why the upper tail is wide even under the strict counterfactual
#
# Under the assumed strict counterfactual we still expect a
# substantial pipeline from sources who became symptomatic just
# before the assumed intervention date. Their offspring chains were
# almost certainly initiated before the cut-off (the fitted
# transmission timing δ is centred near zero days post-source-onset),
# and most of those offspring would still be in the incubation period
# at `obs_date`. The wide upper band reflects two compounding sources
# of uncertainty: the R(t) posterior at the latest knots is wide
# because few cases inform it directly, and the small posterior of
# the dispersion `k` means even sources with no observed offspring
# carry meaningful pipeline mass. See the
# [Limitations page](limitations.md) for the real-time fitting
# caveats.
