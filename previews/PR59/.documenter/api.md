
# API reference {#API-reference}
<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.diagnostics_table-Tuple{Any}' href='#TransmissionLinelist.diagnostics_table-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.diagnostics_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
diagnostics_table(chn)
```


Single-row `DataFrame` summarising sampler diagnostics: maximum R̂, minimum bulk ESS, divergence count, and wall-clock sampling time in seconds. The runtime is read from FlexiChains&#39; per-chain `sampling_time` metadata; under `MCMCThreads` chains run in parallel so the wall clock is approximated by the maximum over chains. Returns `missing` for the runtime if the chain carries no timing metadata.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/plots.jl#L171-L180" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_data-Tuple{Any}' href='#TransmissionLinelist.plot_data-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.plot_data</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_data(ll)
```


Two-panel view of the raw line list: epicurve by ISO week of onset (left) and exposure windows against onset dates (right). Returns a `Makie.Figure`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/plots.jl#L13-L18" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_delta_sense_check-Tuple{Any, Any}' href='#TransmissionLinelist.plot_delta_sense_check-Tuple{Any, Any}'><span class="jlbinding">TransmissionLinelist.plot_delta_sense_check</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_delta_sense_check(chn, data)
```


Sense-check the per-pair posterior of δ against the fitted population `Normal(μ_δ, σ_δ)`. For each sourced pair, take the posterior of `δ_pair = T_inf[secondary] − T_onset[source]` and reduce to its median; then plot the histogram of those per-pair medians with the population density overlaid. Returns a `Makie.Figure`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/plots.jl#L384-L392" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_inc_sense_check-Tuple{Any, Any}' href='#TransmissionLinelist.plot_inc_sense_check-Tuple{Any, Any}'><span class="jlbinding">TransmissionLinelist.plot_inc_sense_check</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_inc_sense_check(chn, data; n_density_draws = 200)
```


Sense-check the per-case posterior of the incubation period against the fitted population `LogNormal(μ_inc, σ_inc)`. For each case, takes the posterior of `inc_i = T_onset[i] − T_inf[i]` and reduces to its median; plots the histogram of those per-case medians with the median PDF (and 95% pointwise ribbon) of the population LogNormal overlaid. Returns a `Makie.Figure`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/plots.jl#L618-L627" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_pair-Tuple{Any}' href='#TransmissionLinelist.plot_pair-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.plot_pair</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_pair(chn; thin = 2)
```


Corner plot of the population scalars `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`, `k` via PairPlots.jl. Returns a Makie `Figure` (requires a Makie backend such as CairoMakie loaded at the call site).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/plots.jl#L193-L199" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_predictive_distributions-Tuple{Any}' href='#TransmissionLinelist.plot_predictive_distributions-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.plot_predictive_distributions</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_predictive_distributions(chn; rng = Random.MersenneTwister(1))
```


Two-by-two panel of the implied population distributions under the posterior for incubation period, transmission timing δ, generation interval, and serial interval. Each panel shows draws from `p(y_new | data) = ∫ p(y_new | θ) p(θ | data) dθ`, i.e. what a new case or transmission pair would look like under the fitted parameters.

This is _not_ a posterior-predictive check against observed data; for that, see `plot_z_ppc`, `plot_delta_sense_check`, and `plot_inc_sense_check`.

Inc and δ panels overlay the parametric density (median PDF with a 95% pointwise ribbon across draws) and a histogram of one predictive sample per draw. GI and SI show the predictive-sample histogram only. Returns a `Makie.Figure`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/plots.jl#L266-L283" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_prior_predictives-Tuple{}' href='#TransmissionLinelist.plot_prior_predictives-Tuple{}'><span class="jlbinding">TransmissionLinelist.plot_prior_predictives</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_prior_predictives(; n = 5000, rng = Random.MersenneTwister(0))
```


Prior-predictive panel: histograms of Inc, δ, and GI/SI drawn from the package&#39;s independent priors on `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`. Returns a `Makie.Figure`.

Three histograms faceted by quantity is the kind of plot AoG was built for: one long-form data frame, `mapping(:value, layout = :panel)`, `visual(Hist)`. Each panel still has its own viewing window so long tails don&#39;t squash the bars; rather than per-facet axis limits, the input is pre-clipped to the window for each panel.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/plots.jl#L667-L679" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_rt-Tuple{Any}' href='#TransmissionLinelist.plot_rt-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.plot_rt</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_rt(chn; n_draws_plot = 100, ymax = 4.0)
```


Spaghetti plot of R(t) over the weekly knots. Each thinned posterior draw is a piecewise-linear trajectory through `(knot_date[b], exp(log_R[b]))`. Knot dates come from `BIN_EDGES` (data.jl). Returns a `Makie.Figure`.

Per-draw spaghetti is built as a long-form `DataFrame` and drawn via AlgebraOfGraphics with `group = :draw`, which is the idiomatic way to spell &quot;one line per draw&quot; once the data is tidy.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/plots.jl#L333-L343" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.sample_fit-Tuple{Any}' href='#TransmissionLinelist.sample_fit-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.sample_fit</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sample_fit(model; samples=1000, chains=4, target_accept=0.95,
           seed=20260508, progress=false)
```


Run NUTS on `model` using the package&#39;s default Enzyme AD backend and `InitFromPrior()` chain initialisation. Returns the FlexiChain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/main.jl#L4-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.summary_table-Tuple{Any}' href='#TransmissionLinelist.summary_table-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.summary_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
summary_table(chn)
```


Posterior summary `DataFrame` for the headline quantities: incubation mean, 95th and 99th percentiles, transmission timing μ_δ / σ_δ, GI / SI mean and SD, and Negative-Binomial dispersion k. Columns: `quantity`, `median`, `lower_95`, `upper_95`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/plots.jl#L123-L130" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.z_ppc_summary-Tuple{Any, Any}' href='#TransmissionLinelist.z_ppc_summary-Tuple{Any, Any}'><span class="jlbinding">TransmissionLinelist.z_ppc_summary</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
z_ppc_summary(chn, d; rng = Random.MersenneTwister(1),
              edges = bin_edges_day(d.t0))
```


Companion to `plot_z_ppc` returning a `DataFrame` of numeric posterior-predictive summaries for three discrete test statistics — `sum(Z)`, `max(Z)`, and `count(Z = 0)`. Replicates `Z_rep` jointly in `(T_inf, log_R, k)` to match `plot_z_ppc`. Columns: `statistic`, `observed`, `rep_median`, `rep_lower_95`, `rep_upper_95`, `p_ppp`, where `p_ppp = 2 · min(P(T_rep ≥ T_obs), P(T_rep ≤ T_obs))` is the two-sided Bayesian posterior-predictive p-value.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/bff0eab80cee6144443bf9c4e27e7eba671f2334/src/plots.jl#L494-L505" target="_blank" rel="noreferrer">source</a></Badge>

</details>

