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
# This page validates the corrections by fitting the same outbreak at four real-time cut-offs and overlaying the resulting R(t) posteriors and population marginals against a counterfactual retrospective and the full closed-out fit.
# It also runs a **delays-only diagnostic** at each cut-off — fitting just the incubation and δ submodels — so that if the full joint fit collapses, we can tell whether the pathology lives in the delay submodels or in the R(t) / `case_model` half of the likelihood.
#
# Real-time-specific caveats are listed in the [Limitations page](limitations.md)
# under the "Real-time fitting caveats" heading.

using TransmissionLinelist
using AlgebraOfGraphics: data, mapping, visual, draw
using DataFrames: DataFrame, groupby, nrow
using DataFramesMeta: @chain, @rtransform, @rsubset, @select, @transform,
                      @combine, @subset
using Dates: Dates, Date, Day
using Random
using Statistics: quantile
using CairoMakie
using CairoMakie: Band, Lines, Hist, VLines
using FlexiChains: FlexiChains
using Turing: DynamicPPL
using Logging: Logging

Random.seed!(20260512)
## Silence the NUTS "Found initial step size" Info logs that would
## otherwise clutter the rendered example output.
Logging.disable_logging(Logging.Info)

# ## Setup
#
# Four cut-offs are used: 31 December 2018 (about three weeks into the
# outbreak), 7 January 2019, 14 January 2019, and 21 January 2019.
# Together they show successive stages of inference as the
# offspring-completeness window expands: at Dec 31 most late-December
# sources have had little time to seed offspring; at Jan 7 the
# completeness factor for those sources has begun to fill in; by
# Jan 14 the late-December chains have largely played out; and by
# Jan 21 the early-January chains have also matured, so R(t) at
# those knots should sharpen toward the retrospective.

obs_dates = [Date("2018-12-31"), Date("2019-01-07"),
    Date("2019-01-14"), Date("2019-01-21")]
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
    fit_dly(d) = sample_fit(delays_only_model(d); seed = seed)
    merge(prep, (; chn_dly_truth = fit_dly(prep.d_truth),
        chn_dly_rt = fit_dly(prep.d_rt)))
end;

# ## Delays-only population posteriors per cut-off
#
# Overlaid marginals of `(μ_inc, σ_inc, μ_δ, σ_δ)` from the delays-only
# fits at each cut-off, one row per `obs_date`. If the corrected
# real-time density tracks the counterfactual retro density on each
# panel, the right-truncation on Inc and δ alone is enough to recover
# the population delays from the truncated observations.

# Shared marginal-overlay helper used by both this delays-only plot
# and the joint-fit version further down.
function plot_marginal_overlay(df; size_kw = (1500, 1200))
    spec = data(df) *
           mapping(:value => "value", color = :fit => "fit",
               row = :obs_date, col = :param) *
           visual(Hist; bins = 30, normalization = :pdf, alpha = 0.4)
    return draw(spec;
        facet = (linkxaxes = :colwise, linkyaxes = :none),
        figure = (; size = size_kw))
end

# Long-form DataFrame of scalar parameter draws from a FlexiChain,
# selecting on `params` via FlexiChains' VarName indexing.
function chain_long(chn, params)
    vns = [DynamicPPL.VarName{p}() for p in params]
    @chain DataFrame(FlexiChains.Long(chn[vns])) begin
        @rtransform :param = String(DynamicPPL.getsym(:param))
        @select :iter :chain :param :value
    end
end

let
    params = [:μ_inc, :σ_inc, :μ_δ, :σ_δ]
    df = @chain delays_fits begin
        map(_) do fit
            t = @transform(chain_long(fit.chn_dly_truth, params),
                :obs_date=string(fit.obs_date),
                :fit="counterfactual retro")
            r = @transform(chain_long(fit.chn_dly_rt, params),
                :obs_date=string(fit.obs_date),
                :fit="corrected real-time")
            vcat(t, r)
        end
        reduce(vcat, _)
    end
    plot_marginal_overlay(df)
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
    fit_jnt(d) = sample_fit(joint_model(d, prep.edges_rt); seed = seed)
    chn_truth = fit_jnt(prep.d_truth)
    chn_rt = fit_jnt(prep.d_rt)
    merge(prep,
        (; m_truth = joint_model(prep.d_truth, prep.edges_rt),
            m_rt = joint_model(prep.d_rt, prep.edges_rt),
            chn_truth, chn_rt,
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
    merge((; obs_date, fit_kind, N_cases = n_cases),
        first(diagnostics_table(chn)))
end

diag_df = let
    full = [diag_row(chn_retro_full, missing,
        "full retro (joint)", nrow(ll))]
    per_cut = mapreduce(vcat, zip(delays_fits, joint_fits)) do (dly, jnt)
        [
            diag_row(dly.chn_dly_truth, dly.obs_date,
                "retro (delays only)", dly.n_truth),
            diag_row(dly.chn_dly_rt, dly.obs_date,
                "realtime (delays only)", dly.n_rt),
            diag_row(jnt.chn_truth, jnt.obs_date,
                "retro (joint)", jnt.n_truth),
            diag_row(jnt.chn_rt, jnt.obs_date,
                "realtime (joint)", jnt.n_rt)]
    end
    DataFrame(vcat(full, per_cut))
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
                    color = :fit, layout = :obs_date) *
                visual(Band; alpha = 0.2)
    line_spec = data(df) *
                mapping(:bin => "Bin index",
                    :med => "R(t) (80% CrI)";
                    color = :fit, layout = :obs_date) *
                visual(Lines; linewidth = 2)
    draw(band_spec + line_spec;
        axis = (; limits = (nothing, (0.0, 4.0))),
        facet = (; linkxaxes = :all, linkyaxes = :all),
        figure = (; size = (1800, 1100)))
end

# ## Population posteriors per cut-off
#
# Overlaid marginals of `(μ_inc, σ_inc, μ_δ, σ_δ, k)` at each cut-off,
# one row per `obs_date`. If the corrected real-time density tracks
# the counterfactual retro density on each panel, the corrections are
# recovering the population from the truncated observations.

# Long-form DataFrame of scalar parameter draws from a `summarise`
# posterior NamedTuple (one Vector{Float64} per param).
function post_long(post, params; obs_date, fit)
    rows = mapreduce(vcat, params) do p
        vals = getproperty(post, p)
        DataFrame(obs_date = string(obs_date), fit = fit,
            param = String(p), value = collect(vals))
    end
    return rows
end

let
    params = [:μ_inc, :σ_inc, :μ_δ, :σ_δ, :k]
    df = @chain joint_fits begin
        map(_) do fit
            vcat(
                post_long(fit.post_truth, params;
                    obs_date = fit.obs_date,
                    fit = "counterfactual retro"),
                post_long(fit.post_rt, params;
                    obs_date = fit.obs_date,
                    fit = "corrected real-time"),
                post_long(post_retro_full, params;
                    obs_date = fit.obs_date,
                    fit = "full retrospective"))
        end
        reduce(vcat, _)
    end
    ## Cap `k` at its 99% quantile (pooled across fits) so the long
    ## right tail doesn't compress the other panels visually.
    k_cap = @chain df begin
        @rsubset :param == "k" && isfinite(:value)
        quantile(_.value, 0.99)
    end
    df_capped = @rsubset(df, :param != "k" || :value <= k_cap)
    plot_marginal_overlay(df_capped; size_kw = (1900, 1500))
end

# ## Controlled-outbreak projection
#
# At each `obs_date`, three counterfactual predictions for the number
# of future symptomatic cases, each conditional on an *assumed*
# intervention rule. None of these scenarios is a factual claim about
# what happened in Epuyén; they describe what the fitted model implies
# *if* transmission had been interrupted in one of three ways.
#
# - **Intervention on 2018-12-31** — transmission is assumed to have
#   stopped on the onset date of the 18th case.
# - **Intervention at the cut-off** — transmission is assumed to have
#   stopped at each panel's `obs_date`. This is
#   [`predict_controlled_outbreak`](@ref)'s default; for the Dec 31
#   panel it coincides with the variant above.
# - **No intervention** — current sources keep transmitting; see
#   [`predict_natural_chain_outbreak`](@ref).
#
# The realised count of cases with onset strictly after each cut-off
# is overlaid as a vertical reference. It is one realisation drawn
# from the true (unknown) process; the predictive distribution is the
# posterior over what could have happened consistent with the fitted
# model and the assumed intervention rule, so we expect spread around
# the realised count rather than agreement on a point.

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

# Per-source breakdown for the intervention-at-cut-off counterfactual
# at the latest `obs_date`: shows which sources drive the projected
# future onset count, one row per observed source.

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
    kinds = [(:strict_samples, "intervention 2018-12-31"),
        (:at_obs_samples, "intervention at cut-off"),
        (:natural_samples, "no intervention")]
    df_hist = @chain controlled begin
        mapreduce(vcat, _) do c
            mapreduce(vcat, kinds) do (field, label)
                DataFrame(panel = string(c.obs_date), kind = label,
                    value = Float64.(getproperty(c, field)))
            end
        end
        groupby(:panel)
        @transform :cap = quantile(:value, 0.99)
        @rsubset :value <= :cap
        @select :panel :kind :value
    end
    df_vline = @chain controlled begin
        DataFrame(panel = string.(getfield.(_, :obs_date)),
            kind = "realised",
            value = Float64.(getfield.(_, :actual_future)))
    end
    hist_spec = data(df_hist) *
                mapping(:value => "Future cases";
                    color = :kind, col = :panel) *
                visual(Hist; bins = 30, normalization = :pdf, alpha = 0.4)
    vline_spec = data(df_vline) *
                 mapping(:value; color = :kind, col = :panel) *
                 visual(VLines; linewidth = 3)
    draw(hist_spec + vline_spec;
        facet = (linkxaxes = :none, linkyaxes = :none),
        figure = (; size = (1900, 400)))
end

# ### Why the Jan 7 predictive is wider than Dec 31, and how Jan 14 tightens
#
# More data does not always mean a tighter future-cases prediction.
# At Jan 7 the line list adds secondary cases with onsets in the last
# week before the cut-off, so for those sources the chain is far from
# complete and most of their expected offspring still lies in the
# future.
# At Dec 31 the late-onset sources are fewer and most of their
# expected offspring chains have already completed, so each source
# contributes less to the future window.
# By Jan 14 the late-December sources have had two further weeks for
# their offspring to surface, so the completeness factor pins their
# per-source rate more tightly and the predictive distribution
# narrows even though the late-cut-off sources themselves remain
# pipeline-heavy.
# By Jan 21 the same is true of the early-January cohort, so the
# predictive tightens further toward the realised count.
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
#
# !!! note "Per-source intervention rules"
#
#     The scenarios above set a single intervention date (or none)
#     that applies to every source. `intervention_time` also accepts
#     a `Vector{Date}` of length `N` — one cut-off per observed
#     source — for cases where isolation policies differ across
#     individuals. An aggressive contact-tracing scenario in which
#     each case is isolated on or shortly after their own symptom
#     onset could be encoded by passing a per-source vector derived
#     from `ll_rt.onset_date`. We do not have direct evidence of
#     per-case isolation in the Epuyén outbreak so we do not run
#     this scenario here, but the API supports it.
