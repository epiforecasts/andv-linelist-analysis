# # Intervention-aware R(t) walkthrough
#
# The base `joint_model` fits R(t) as an unconstrained weekly random walk.
# Outbreaks with a known intervention may invite stronger structure: a
# break point, or a random walk plus a one-sided shock that can only
# drop R after the intervention date.
# This page compares three intervention-aware variants on the same
# real-time cut-off — the 18-case wave at `Date("2018-12-31")` —
# holding everything else constant.
#
# All three scenarios are fit-side assumptions about R(t). They are not
# claims about what happened in Epuyén; they describe what the data
# imply *if* we constrain R(t) to behave a particular way around the
# assumed intervention date.

using TransmissionLinelist
using AlgebraOfGraphics: data, mapping, visual, draw
using DataFrames: DataFrame, groupby, nrow
using DataFramesMeta: @chain, @rtransform, @rsubset, @select, @transform,
                      @combine, @subset
using Dates: Dates, Date, Day
using Distributions: Distributions, Normal, LogNormal, NegativeBinomial,
                     truncated, logpdf, cdf
using Random
using Statistics: quantile, mean, var
using CairoMakie
using CairoMakie: Band, Lines, Hist
using FlexiChains: FlexiChains
using Turing: Turing, @model, to_submodel, DynamicPPL
using Logging: Logging

Random.seed!(20260512)
## Silence the NUTS "Found initial step size" Info logs.
Logging.disable_logging(Logging.Info)

# ## Setup
#
# Load the closed-out line list and build the model inputs and weekly
# knot grid. The intervention timing is a model-side argument, not a
# data cut-off — the fits see the full outbreak.

intervention_date = Date("2018-12-31")
ll = load_linelist();
t0_ref = minimum(ll.onset_date) - Day(60)
d = build_data(ll; t0 = t0_ref)
edges = prepare_rt_edges(t0_ref)
intervention_day = Float64(Dates.value(intervention_date - t0_ref))
seed = 20260512

# ## Submodel variants
#
# Two new R(t) submodels are defined inline for this page. Both return
# `(; log_R)` so the joint model's `case_model` consumes them
# unchanged. `log_R` is registered on the chain via `:=` so downstream
# summaries (`vector_chain(chn, :log_R)`, `rt_band`, …) keep working.

# ### Two-level step (intervention-only)
#
# A single break point at `intervention_day`: one log-R before, one
# after. Each knot's log-R is set to one of the two scalars depending
# on which side of the break point it falls.

@model function step_rt_submodel(edges, intervention_day::Real;
        pre_prior = Normal(log(1.5), 1.0),
        post_prior = Normal(log(1.0), 1.0))
    log_R_pre ~ pre_prior
    log_R_post ~ post_prior
    log_R := [e < intervention_day ? log_R_pre : log_R_post for e in edges]
    return (; log_R)
end

# ### Random walk plus one-sided post-intervention shock
#
# Standard random walk plus a multiplicative `exp(shock)` applied at
# every knot on or after `intervention_day`. The shock prior is
# truncated above at zero so the model can only express a drop in R
# after the intervention.

@model function rw_plus_shock_rt_submodel(edges, intervention_day::Real;
        init_prior = Normal(log(1.5), 1.0),
        sigma_prior = truncated(Normal(0.0, 0.2); lower = 0),
        shock_prior = truncated(Normal(0.0, 0.5); upper = 0))
    n_knots = length(edges)
    σ_rw ~ sigma_prior
    log_R_init ~ init_prior
    shock ~ shock_prior
    T = typeof(log_R_init)
    ε ~ Turing.filldist(Normal(zero(T), one(T)), n_knots - 1)
    base = vcat(log_R_init, log_R_init .+ accumulate(+, σ_rw .* ε))
    log_R := [base[i] + (edges[i] >= intervention_day ? shock : zero(T))
              for i in 1:n_knots]
    return (; log_R)
end

# ## Fits
#
# Three joint fits on the same real-time view, identical seed,
# differing only in the R(t) submodel.

m_rw = joint_model(d, edges)
m_step = joint_model(d, edges;
    rt = step_rt_submodel(edges, intervention_day))
m_shock = joint_model(d, edges;
    rt = rw_plus_shock_rt_submodel(edges, intervention_day))

chn_rw = sample_fit(m_rw; seed = seed)
chn_step = sample_fit(m_step; seed = seed)
chn_shock = sample_fit(m_shock; seed = seed)

post_rw = summarise(chn_rw)
post_step = summarise(chn_step)
post_shock = summarise(chn_shock)

# ## Sampler health
#
# R̂ near 1 and zero divergences are the targets across all three fits.
# Step and RW + shock have one extra scalar each (`log_R_post`,
# `shock`) so a small uptick in `runtime_seconds` is expected.

diag_df = let
    rows = []
    for (name, chn) in (("RW (current)", chn_rw),
        ("step (intervention only)", chn_step),
        ("RW + post-shock", chn_shock))
        push!(rows, merge((; scenario = name),
            first(diagnostics_table(chn))))
    end
    DataFrame(rows)
end

# ## R(t) trajectories overlaid
#
# Posterior median R(t) with 80% CrI ribbons for all three scenarios on
# a single panel. The step fit forces a flat-flat shape; the RW + shock
# fit lets the random walk move freely up to `intervention_day` and
# only allows a downward jump at the cut-off; the plain RW fit imposes
# no intervention structure.

let
    band_rows = DataFrame[]
    for (name, post) in (("RW (current)", post_rw),
        ("step (intervention only)", post_step),
        ("RW + post-shock", post_shock))
        tbl = rt_band(post)
        tbl.scenario .= name
        push!(band_rows, tbl)
    end
    df = reduce(vcat, band_rows)
    band_spec = data(df) *
                mapping(:bin => "Bin index",
                    :lo => "R(t) (80% CrI)", :hi;
                    color = :scenario) *
                visual(Band; alpha = 0.2)
    line_spec = data(df) *
                mapping(:bin => "Bin index",
                    :med => "R(t) (80% CrI)";
                    color = :scenario) *
                visual(Lines; linewidth = 2)
    draw(band_spec + line_spec;
        axis = (; limits = (nothing, (0.0, 4.0))),
        figure = (; size = (900, 500)))
end

# ## Population marginal posteriors
#
# Overlaid marginals of `(μ_inc, σ_inc, μ_δ, σ_δ, k)` across the three
# scenarios. If the population delays move much when only the R(t)
# half of the likelihood changes, that hints the data carries little
# information on the delays at this cut-off; otherwise they should
# overlap.

function plot_marginal_overlay(df; size_kw = (1600, 400))
    spec = data(df) *
           mapping(:value => "value", color = :scenario => "scenario",
               col = :param) *
           visual(Hist; bins = 30, normalization = :pdf, alpha = 0.4)
    return draw(spec;
        facet = (linkxaxes = :colwise, linkyaxes = :none),
        figure = (; size = size_kw))
end

function post_long(post, params; scenario)
    rows = mapreduce(vcat, params) do p
        vals = getproperty(post, p)
        DataFrame(scenario = scenario, param = String(p),
            value = collect(vals))
    end
    return rows
end

let
    params = [:μ_inc, :σ_inc, :μ_δ, :σ_δ, :k]
    df = vcat(
        post_long(post_rw, params; scenario = "RW (current)"),
        post_long(post_step, params;
            scenario = "step (intervention only)"),
        post_long(post_shock, params; scenario = "RW + post-shock"))
    k_cap = @chain df begin
        @rsubset :param == "k" && isfinite(:value)
        quantile(_.value, 0.99)
    end
    df_capped = @rsubset(df, :param != "k" || :value <= k_cap)
    plot_marginal_overlay(df_capped)
end

# ## R-pre vs R-post in the intervention-bearing scenarios
#
# For the step model, R-pre and R-post are the two sampled scalars
# directly. For the RW + post-shock model, R-pre is the exponentiated
# log-R at the last pre-intervention knot and R-post is R-pre times
# `exp(shock)`; both pulled from the per-draw `:log_R` vector to stay
# consistent with the `case_model` likelihood.

function effect_size_rows(post; scenario, intervention_idx)
    pre = exp.(post.log_R_chain[intervention_idx - 1])
    post_R = exp.(post.log_R_chain[intervention_idx])
    return DataFrame(
        scenario = scenario,
        median_pre = quantile(pre, 0.5),
        lo_pre = quantile(pre, 0.1),
        hi_pre = quantile(pre, 0.9),
        median_post = quantile(post_R, 0.5),
        lo_post = quantile(post_R, 0.1),
        hi_post = quantile(post_R, 0.9),
        p_drop = mean(post_R .< pre))
end

## Identify the first knot at or after `intervention_day`. With the
## real-time grid the cut-off is the last edge, so this is the final
## knot index.
intervention_idx = findfirst(e -> e >= intervention_day, edges)

effect_df = vcat(
    effect_size_rows(post_step;
        scenario = "step (intervention only)",
        intervention_idx = intervention_idx),
    effect_size_rows(post_shock;
        scenario = "RW + post-shock",
        intervention_idx = intervention_idx))

# ## WAIC comparison on the offspring likelihood
#
# Cheap model comparison: pointwise log-likelihood of the observed
# offspring count `Z_i` under each posterior draw, then WAIC =
# `lppd - p_waic` summed across cases. Higher (less negative) is
# better. This conditions on the latent `T_onset` and `log_R` per
# draw, so it is a conditional WAIC over the case-model half of the
# likelihood — the cleanest cross-scenario comparison given that the
# R(t) submodel is the only thing changing.

function pointwise_z_loglik(chn, d, edges)
    log_R = vector_chain(chn, :log_R)
    T_onset = vector_chain(chn, :T_onset)
    k_draws = vec(collect(chn[:k]))
    n_draws = length(k_draws)
    N = d.N
    ll = Matrix{Float64}(undef, n_draws, N)
    p_per_source = if d.obs_time !== nothing
        ## Reproduce the offspring-completeness thinning so the
        ## per-case rate matches the one inside `case_model`.
        μ_inc = vec(collect(chn[:μ_inc]));
        σ_inc = vec(collect(chn[:σ_inc]))
        μ_δ = vec(collect(chn[:μ_δ]));
        σ_δ = vec(collect(chn[:σ_δ]))
        function compute_p(dr, i)
            inc = LogNormal(μ_inc[dr], σ_inc[dr])
            δd = Normal(μ_δ[dr], σ_δ[dr])
            cdf(ConvolvedDelays(inc, δd),
                d.obs_time[i] - T_onset[i][dr])
        end
        compute_p
    else
        (dr, i) -> 1.0
    end
    for dr in 1:n_draws
        logR_d = [log_R[b][dr] for b in eachindex(log_R)]
        for i in 1:N
            R_i = exp(log_R_at(T_onset[i][dr], edges, logR_d))
            p_i = p_per_source(dr, i)
            R_eff = R_i * p_i
            kd = k_draws[dr]
            p = max(kd / (kd + R_eff), eps(typeof(kd)))
            ll[dr, i] = logpdf(NegativeBinomial(kd, p), d.Zobs[i])
        end
    end
    return ll
end

function waic(ll)
    ## Standard WAIC: lppd − p_waic, summed across observations.
    lppd = sum(log.(mean(exp.(ll); dims = 1)))
    p_waic = sum(var(ll; dims = 1))
    return (; waic = lppd - p_waic, lppd, p_waic)
end

waic_df = let
    rows = []
    for (name, chn) in (("RW (current)", chn_rw),
        ("step (intervention only)", chn_step),
        ("RW + post-shock", chn_shock))
        ll = pointwise_z_loglik(chn, d, edges)
        w = waic(ll)
        push!(rows, (; scenario = name, w...))
    end
    DataFrame(rows)
end

# ## Synthesis
#
# The plain RW fit and the RW + post-shock fit both put substantial
# posterior mass on R(t) below 1 in the final week of the cut-off
# without requiring an explicit break point, so the data alone are
# already nudging late R downward.
# The step fit forces a discrete break, which buys identifiability of
# R-pre and R-post but cannot represent the gradual rise-and-fall that
# the RW fit picks out.
# The RW + post-shock fit nests the plain RW (`shock = 0` is at the
# upper boundary of the shock prior) so its `log_R` posterior can only
# move down at the break: if WAIC barely improves over the plain RW,
# the data are not asking for the extra shock parameter.
# The R-pre vs R-post tables in `effect_df` quantify the drop the data
# attribute to the assumed intervention under each parameterisation;
# the headline `p_drop` summarises the posterior probability that
# R-post is strictly below R-pre.
#
# A fourth scenario varying the right-truncation correction on the
# delay distributions is discussed in
# [issue #48](https://github.com/epiforecasts/andv-linelist-analysis/issues/48).
