# # Prior sensitivity
#
# The joint model ships with a default set of priors: `μ_inc`,
# `σ_inc` on the incubation lognormal; `μ_δ`, `σ_δ` on the
# transmission timing normal; `1/√k` on the offspring NB dispersion;
# and `σ_rw` on the weekly random walk for log R(t). This page sweeps
# one prior at a time and overlays the resulting marginal posteriors
# so the influence of each prior on the headline parameters is
# visible at a glance.
#
# Each scenario is a single `joint_model(...)` call with one submodel
# kwarg overridden, fitted with `sample_fit` under a shared seed.
# Posterior summaries are read directly from the fitted chains; no
# numbers are hard-coded.

using TransmissionLinelist
using AlgebraOfGraphics: data, mapping, visual, draw
using DataFrames: DataFrame, nrow
using DataFramesMeta: @chain, @rtransform, @rsubset, @select, @transform
using Dates: Date, Day
using Distributions: Normal, truncated
using Random
using Statistics: quantile
using CairoMakie
using CairoMakie: Band, Lines
using FlexiChains: FlexiChains
using Turing: DynamicPPL
using Logging: Logging

Random.seed!(20260514)
## Silence the NUTS "Found initial step size" Info logs that would
## otherwise clutter the rendered example output.
Logging.disable_logging(Logging.Info)

# ## Setup
#
# A single retrospective line list provides the data. R(t) knots use
# the package default. All scenarios share the same data and seed; the
# only difference between them is the prior on one submodel.

ll = load_linelist();
t0_ref = minimum(ll.onset_date) - Day(60)
d_retro = build_data(ll; t0 = t0_ref)
edges_ref = prepare_rt_edges(t0_ref)
seed = 20260514

# ## Scenarios
#
# Five scenarios, all fit to the closed-out outbreak:
#
# - **baseline** — package defaults; `σ_rw ~ N⁺(0, 0.2)`,
#   `μ_inc ~ N(3.0, 0.5)`, `σ_δ ~ N⁺(0, 1.0)`, `1/√k ~ N⁺(0, 1.0)`.
# - **wider σ_rw** — the old default `σ_rw ~ N⁺(0, 0.5)`, allowing
#   week-on-week R(t) jumps that are 2.5× larger a priori.
# - **tighter Inc μ** — `μ_inc ~ N(3.0, 0.1)` instead of `(3.0, 0.5)`;
#   the population log-mean incubation is pinned closer to 3.0.
# - **k prior shifted** — `1/√k ~ N⁺(0.5, 0.5)`, pulling the prior on
#   `k` toward smaller values (more overdispersion).
# - **δ wider** — `σ_δ ~ N⁺(0, 3.0)`, widening the prior on the
#   per-pair transmission-timing SD.

scenarios = [
    (name = "baseline", model = joint_model(d_retro, edges_ref)),
    (name = "wider σ_rw",
        model = joint_model(d_retro, edges_ref;
            rt = random_walk_rt_model(length(edges_ref);
                sigma_prior = truncated(Normal(0.0, 0.5); lower = 0)))),
    (name = "tighter Inc μ",
        model = joint_model(d_retro, edges_ref;
            incubation = incubation_model(;
                μ_prior = Normal(3.0, 0.1)))),
    (name = "k prior shifted",
        model = joint_model(d_retro, edges_ref;
            dispersion = nb_dispersion_model(
                truncated(Normal(0.5, 0.5); lower = 0)))),
    (name = "δ wider",
        model = joint_model(d_retro, edges_ref;
            transmission = transmission_delta_model(;
                σ_prior = truncated(Normal(0.0, 3.0); lower = 0))))
]

# ## Fitting
#
# All scenarios share the same NUTS settings (default warmup,
# 1000 samples × 4 chains) and the same seed, so any difference in the
# overlaid marginals is attributable to the prior swap.

fits = map(scenarios) do s
    chn = sample_fit(s.model; seed = seed)
    (; s.name, chn, post = summarise(chn))
end;

# ## Sampler health
#
# One row per scenario. R̂ near 1 and zero divergences are the targets;
# if a wider prior pushes the sampler into a bad geometry it shows up
# here before the marginals are even plotted.

diag_df = DataFrame(
    map(fits) do f
    merge((; scenario = f.name), first(diagnostics_table(f.chn)))
end)

# ## Marginal posteriors
#
# Overlaid marginals of `(μ_inc, σ_inc, μ_δ, σ_δ, k)` across scenarios.
# One row per scenario, one column per parameter; the baseline density
# is repeated as a fixed reference inside every row so each row is a
# direct overlay of one variant against the baseline.

# Long-form DataFrame of scalar parameter draws from one `summarise`
# posterior NamedTuple; mirrors the helper used in the real-time page.
function post_long(post, params; scenario, fit)
    rows = mapreduce(vcat, params) do p
        vals = getproperty(post, p)
        DataFrame(scenario = scenario, fit = fit,
            param = String(p), value = collect(vals))
    end
    return rows
end

let
    params = [:μ_inc, :σ_inc, :μ_δ, :σ_δ, :k]
    baseline = fits[1]
    df = @chain fits begin
        map(_) do f
            vcat(
                post_long(f.post, params;
                    scenario = f.name, fit = f.name),
                post_long(baseline.post, params;
                    scenario = f.name, fit = "baseline (reference)"))
        end
        reduce(vcat, _)
    end
    ## Cap `k` at its 99% quantile (pooled across scenarios) so the
    ## long right tail under the shifted-k prior doesn't compress the
    ## other panels visually.
    k_cap = @chain df begin
        @rsubset :param == "k" && isfinite(:value)
        quantile(_.value, 0.99)
    end
    df_capped = @rsubset(df, :param != "k" || :value <= k_cap)
    plot_marginal_overlay(df_capped;
        row_col = :scenario, col_col = :param,
        size_kw = (1900, 1500))
end

# ## R(t) under each scenario
#
# Posterior medians with 80% CrI ribbons over the weekly knots, one
# colour per scenario. The widely shared knots make it easy to read
# off where the random-walk prior matters most: at the start of the
# outbreak (few cases) and the tail (few completed offspring chains).

let
    band_rows = DataFrame[]
    for f in fits
        tbl = rt_band(f.post)
        tbl.scenario .= f.name
        push!(band_rows, tbl)
    end
    df = reduce(vcat, band_rows)
    band_spec = data(df) *
                mapping(:bin => "Bin index",
                    :lo => "R(t) (80% CrI)", :hi;
                    color = :scenario) *
                visual(Band; alpha = 0.15)
    line_spec = data(df) *
                mapping(:bin => "Bin index",
                    :med => "R(t) (80% CrI)";
                    color = :scenario) *
                visual(Lines; linewidth = 2)
    draw(band_spec + line_spec;
        axis = (; limits = (nothing, (0.0, 4.0))),
        figure = (; size = (1100, 600)))
end

# ## Where the priors actually bite
#
# Reading across the panels:
#
# - **`σ_rw` prior.** Widening to `N⁺(0, 0.5)` mostly inflates the R(t)
#   ribbon at the edges of the knot grid. Headline scalars
#   (`μ_inc`, `μ_δ`, `k`) shift very little, since the random walk
#   prior is essentially independent of the delay submodels in the
#   likelihood.
# - **Inc μ prior.** A tight `N(3.0, 0.1)` on `μ_inc` collapses the
#   `μ_inc` marginal onto its prior mean; `σ_inc` then absorbs the
#   slack and shifts upward to keep the fitted lognormal close to the
#   data. δ and k are largely unaffected.
# - **k prior shift.** Pushing the prior on `1/√k` toward larger
#   values pulls posterior `k` toward smaller values (more dispersion)
#   and widens its right tail. Delay marginals are insensitive.
# - **δ wider.** Allowing a larger `σ_δ` mainly fattens the `σ_δ`
#   posterior itself; `μ_δ` and the incubation parameters barely
#   move, because the per-pair likelihood already pins δ tightly.
#
# The takeaway is that the joint model's headline incubation and
# δ summaries are robust to the exact shape of their submodel priors
# within the ranges swept here. The two priors that matter most for
# downstream R(t) inference are the random-walk SD (which sets the
# smoothness of the R(t) curve) and the dispersion `k` prior (which
# governs how much pipeline mass each observed source can contribute
# without observed offspring). See the
# [Real-time monitoring page](realtime.md) for how these defaults
# behave under real-time cut-offs.
