# # Analysis walkthrough
#
# The Epuyén 2018–19 outbreak in north-west Patagonia was the first cluster where person-to-person Andes hantavirus transmission was documented at scale.
# The line list bundled with this package is hand-encoded from Table S2 of [Martínez et al. 2020](https://doi.org/10.1056/NEJMoa2009040) — 34 cases with exposure windows, symptom-onset dates, and attributed source cases.
#
# This page fits the joint model in `Hantavirus.jl` to that line list and renders the headline outputs.
# Four quantities are estimated together: the incubation period, the transmission timing of each secondary relative to its source's symptom onset (δ), a weekly time-varying reproduction number R(t), and the offspring dispersion `k` of a Negative-Binomial.
# Exposure and onset dates are interval-censored.
# The model handles that by giving each case a continuous latent infection time and a continuous latent onset time, each sampled within its recorded window.
# Generation interval and serial interval are derived in post-processing as the transmission timing plus an incubation period (the source's for GI, the secondary's for SI).
# Fitting all four jointly propagates uncertainty between them that a delay-then-R(t) pipeline would lose.
#
# Priors, the data-augmentation construction, and per-pair GI > 0 constraint are detailed on the [Model](model.md) page.
# Caveats around exposure encoding, late R(t) bins, and right-truncation are on the [Limitations](limitations.md) page.

using Hantavirus
using Chain
using DataFrames
using DataFramesMeta
using Distributions
using FlexiChains
using Plots
using Printf
using Random
using Turing
using ADTypes
using Enzyme
using CairoMakie: CairoMakie
using PairPlots

Random.seed!(20260508)

# ## Load the line list
#
# `load_linelist` parses the bundled CSV and drops the `_alt` sensitivity rows.
# `build_data` re-encodes exposure / onset windows as day offsets from `t0`, 60 days before the first onset.

ll = load_linelist()
d  = build_data(ll)

@chain ll begin
    @select(:patient_id, :exposure_lower, :exposure_upper,
            :onset_date, :source_case, :Z)
    first(8)
end

# ## What the data looks like

plot_data(ll)

# ## Model
#
# The model is `joint_model(d, edges)`.
# Priors are at the top, the random walk on `log R(t)` is non-centred, and the per-pair constraint `T_inf[secondary] > T_inf[source]` rejects trajectories with a non-positive generation interval.

edges = bin_edges_day(d.t0)
model = joint_model(d, edges)

# Source rendered inline so the page never drifts from the implementation.
# The fenced block below is generated from `src/model.jl` at doc-build time via `CodeTracking.@code_string`.

#md # ```@eval
#md # using Markdown
#md # Markdown.MD(Markdown.Code("julia",
#md #     read(joinpath(@__DIR__, "..", "examples", "joint_model_source.jl"), String)))
#md # ```

# ## Fitting
#
# Gradients via Enzyme reverse-mode AD — fastest backend for this model under Turing 0.45.
# Chains are initialised from the prior because the non-centred random walk is sensitive to default zero-initialisation.
# Budget matches the package default: 1000 post-warmup draws across 4 chains, `target_accept = 0.95`.

adtype = AutoEnzyme(; mode = Enzyme.set_runtime_activity(Enzyme.Reverse))

chn = sample(
    model,
    NUTS(0.95; adtype),
    MCMCThreads(),
    1000,
    4;
    initial_params = fill(Turing.DynamicPPL.InitFromPrior(), 4),
    progress = false,
)

diagnostics_table(chn)

# ## Key outputs

summary_table(chn)

# ## R(t) over weekly bins
#
# Spaghetti of thinned posterior draws over the weekly bins; reverts to the prior in late-January bins where cases are thin (see Limitations).

post = mktemp() do _path, io
    redirect_stdout(io) do
        summarise(chn)
    end
end
Hantavirus.plot_rt(post, "rt_spaghetti.png"; n_draws_plot = 100)
nothing #hide

# ![R(t) spaghetti](rt_spaghetti.png)

# ## Pair plot of population parameters
#
# Corner plot for `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`, `k`.

plot_pair(chn)

# ## Posterior-predictive delay distributions
#
# Each panel shows two things.
# The blue band is the posterior over the parametric density: median PDF with a 95% pointwise ribbon across draws.
# The orange histogram is one predictive realisation per draw, pooled.
# GI / SI use a moment-matched Normal for the ribbon and exact `δ + Inc` samples for the histogram.

plot_posterior_predictive(chn)
