# # Analysis walkthrough
#
# The Epuyén 2018–19 outbreak in north-west Patagonia was the first cluster where person-to-person Andes hantavirus transmission was documented at scale.
# The line list bundled with this package is hand-encoded from Table S2 of [Martínez et al. 2020](https://doi.org/10.1056/NEJMoa2009040) — 34 cases with exposure windows, symptom-onset dates, and attributed source cases.
#
# This page fits the joint model in `TransmissionLinelist.jl` to that line list and renders the headline outputs.
# Four quantities are estimated together: the incubation period, the transmission timing of each secondary relative to its source's symptom onset (δ), a weekly time-varying reproduction number R(t), and the offspring dispersion `k` of a Negative-Binomial.
# Exposure and onset dates are interval-censored.
# The model handles that by giving each case a continuous latent infection time and a continuous latent onset time, each sampled within its recorded window.
# Generation interval and serial interval are derived in post-processing as the transmission timing plus an incubation period (the source's for GI, the secondary's for SI).
# Fitting all four jointly propagates uncertainty between them that a delay-then-R(t) pipeline would lose.
#
# Priors, the data-augmentation construction, and per-pair GI > 0 constraint are detailed on the [Model](model.md) page.
# Caveats around exposure encoding, late R(t) bins, and right-truncation are on the [Limitations](limitations.md) page.

using TransmissionLinelist
using Chain
using DataFrames
using DataFramesMeta
using FlexiChains
using Printf
using Random
using CairoMakie
using AlgebraOfGraphics
using PairPlots

Random.seed!(20260508)

# ## Load the line list
#
# `load_linelist` parses the bundled CSV and drops the `_alt` sensitivity rows.
# `build_data` re-encodes exposure / onset windows as day offsets from `t0` (60 days before the first onset); `bin_edges_day` returns the weekly R(t) knot dates as day offsets.

ll = load_linelist()
d = build_data(ll)
edges = bin_edges_day(d.t0)
model = joint_model(d, edges)

@chain ll begin
    @select(:patient_id, :exposure_lower, :exposure_upper,
        :onset_date, :source_case, :Z)
    first(8)
end

# ## What the data looks like

plot_data(ll)

# ## Model

#md # ```@eval
#md # using Markdown
#md # Markdown.MD(Markdown.Code("julia",
#md #     read(joinpath(@__DIR__, "..", "examples", "joint_model_source.jl"), String)))
#md # ```

# ## Prior predictives
#
# Implied prior distributions for the incubation period, transmission timing δ, and the derived generation / serial interval before any data are seen.

plot_prior_predictives()

# ## Fitting
#
# `sample_fit` wraps the package's default NUTS configuration: Enzyme reverse-mode AD, chains initialised from the prior, 1000 post-warmup draws across 4 chains, `target_accept = 0.95`.

chn = sample_fit(model)

# ## Diagnostics
#
# Maximum R̂, minimum bulk ESS, divergence count, and wall-clock sampling time (seconds, approximated by the slowest chain under `MCMCThreads`).

diagnostics_table(chn)

# ## Key outputs

summary_table(chn)

# ## Observed offspring counts
#
# Each case is plotted at its observed offspring count `Z` (with small vertical jitter to break ties) against time.
# Thin segments are posterior draws of the latent infection time `T_inf` joined to the latent onset time `T_onset`; filled dots are the posterior medians of `T_inf` and hollow dots the medians of `T_onset`.
# Index and sourced cases are coloured separately.
# This is a direct view of what the model has to explain before R(t) enters the picture; the posterior-predictive comparison against modelled offspring counts comes later under [`plot_z_ppc`](@ref).

plot_z_dumbbell(chn, d)

# ## R(t) over weekly knots
#
# Spaghetti of thinned posterior draws through the weekly knots (linearly interpolated); reverts to the prior in late-January knots where cases are thin (see Limitations).

plot_rt(chn)

# ## Pair plot of population parameters
#
# Corner plot for `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`, `k`.

plot_pair(chn)

# ## Predictive distributions for the delays
#
# Implied population distributions for the incubation period, transmission timing δ, and the derived generation and serial intervals under the fitted posterior — i.e. what a new case or transmission pair would look like.
# This is not a check against the observed data; for that see the sense-check and PPC panels below.
# Inc and δ panels show the posterior over the parametric density (median PDF with a 95% pointwise ribbon across draws) overlaid with one predictive realisation per draw; GI and SI show the predictive-sample histogram only.

plot_predictive_distributions(chn)

# ## δ sense check
#
# Compare the per-pair posterior medians of δ to the fitted population `Normal(μ_δ, σ_δ)`.

plot_delta_sense_check(chn, d)

# ## Incubation-period sense check
#
# Compare the per-case posterior medians of `T_onset[i] − T_inf[i]` to the fitted population `LogNormal(μ_inc, σ_inc)`.

plot_inc_sense_check(chn, d)

# ## Offspring posterior-predictive check
#
# Joint-draw posterior-predictive check.
# For each posterior draw, replicate `Z_rep[i] ~ NegativeBinomial(k, k/(k+R_i))` per case, with `R_i = exp(log_R_at(T_inf[i], edges, log_R))` evaluated at the same draw's `T_inf[i]`, `log_R`, and `k` (matching the model's likelihood, clamp and all).
# The left panel compares frequencies of each `Z` value against the observed line list.
# The right column has three stacked subpanels — one per discrete test statistic (`sum(Z)`, `max(Z)`, `count(Z = 0)`) — each showing the histogram of the replicated statistic with the observed value as a dashed vertical rule.

plot_z_ppc(chn, d)

# Numeric values for each test statistic — observed, replicated median + 95% CrI, and the two-sided Bayesian posterior-predictive p-value `2 · min(P(T_rep ≥ T_obs), P(T_rep ≤ T_obs))`.

z_ppc_summary(chn, d)
