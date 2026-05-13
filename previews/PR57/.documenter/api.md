
# API reference {#API-reference}
<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.bin_edges_day-Tuple{Any}' href='#TransmissionLinelist.bin_edges_day-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.bin_edges_day</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
bin_edges_day(t0) -> Vector{Float64}

```


Return the weekly R(t) knot dates expressed as days relative to `t0`.

The knots span the outbreak in weekly steps; combined with [`log_R_at`](/api#TransmissionLinelist.log_R_at-Tuple{Real,%20AbstractVector{<:Real},%20Any}) this defines the piecewise-linear log R(t) trajectory used by [`joint_model`](/api#TransmissionLinelist.joint_model-Tuple{Any,%20Any}).

**Arguments**
- `t0`: the model&#39;s time origin (the `t0` field of the tuple returned by [`build_data`](/api#TransmissionLinelist.build_data-Tuple{Any})).
  

**Returns**

A `Vector{Float64}` of length `length(BIN_EDGES)` giving the knot positions in days.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/data.jl#L105" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.bin_labels-Tuple{}' href='#TransmissionLinelist.bin_labels-Tuple{}'><span class="jlbinding">TransmissionLinelist.bin_labels</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
bin_labels() -> Vector{String}

```


Return string labels for the weekly R(t) knots.

One entry per `log_R` element. Used to label plots and posterior summaries produced by [`summarise`](/api#TransmissionLinelist.summarise-Tuple{Any}) and [`plot_rt`](/api#TransmissionLinelist.plot_rt-Tuple{Any}).

**Returns**

A `Vector{String}` of ISO-format knot dates.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/data.jl#L174" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.build_data-Tuple{Any}' href='#TransmissionLinelist.build_data-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.build_data</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
build_data(
    ll
) -> NamedTuple{(:t0, :onset_lo_day, :onset_hi_day, :exp_lo_day, :exp_hi_day, :source_idx, :Zobs, :N), <:NTuple{8, Any}}

```


Build the model input from a line-list `DataFrame`.

Anchors all times in days relative to `t0 = minimum(onset_date) - 60 d` and encodes interval-censored onset and exposure windows as `[lo, hi)` pairs. Resolves the `source_case` column to integer indices into the line list (`0` denotes a zoonotic index case with no human source) and reads observed offspring counts from the `Z` column.

**Arguments**
- `ll`: a `DataFrame` returned by [`load_linelist`](/api#TransmissionLinelist.load_linelist).
  

**Returns**

A named tuple `(t0, onset_lo_day, onset_hi_day, exp_lo_day, exp_hi_day, source_idx, Zobs, N)` ready to pass to [`joint_model`](/api#TransmissionLinelist.joint_model-Tuple{Any,%20Any}).

**Examples**

```julia
ll = load_linelist()
d  = build_data(ll)
d.N
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/data.jl#L61" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.diagnostics-Tuple{Any}' href='#TransmissionLinelist.diagnostics-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.diagnostics</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
diagnostics(
    chn
) -> NamedTuple{(:rhat, :ess, :ndiv), <:Tuple{Float64, Float64, Any}}

```


Return convergence diagnostics for `chn`: `(; rhat, ess, ndiv)` — the maximum `R̂` across scalar parameter entries, the minimum bulk ESS, and the divergent transition count.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/postprocess.jl#L32" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.diagnostics_table-Tuple{Any}' href='#TransmissionLinelist.diagnostics_table-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.diagnostics_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
diagnostics_table(chn)
```


Single-row `DataFrame` summarising sampler diagnostics: maximum R̂, minimum bulk ESS, divergence count, and wall-clock sampling time in seconds. The runtime is read from FlexiChains&#39; per-chain `sampling_time` metadata; under `MCMCThreads` chains run in parallel so the wall clock is approximated by the maximum over chains. Returns `missing` for the runtime if the chain carries no timing metadata.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/plots.jl#L171-L180" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.joint_model-Tuple{Any, Any}' href='#TransmissionLinelist.joint_model-Tuple{Any, Any}'><span class="jlbinding">TransmissionLinelist.joint_model</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
joint_model(d, edges) -> Any

```


Joint Turing model for the Epuyén ANDV outbreak.

Estimates the incubation period (LogNormal), the per-pair transmission timing `δ` relative to source onset (Normal), a weekly piecewise-linear log-R(t) random walk, and the offspring dispersion `k`, from interval- censored exposure and onset windows. Each case has continuous latent infection and onset times; the positive generation-interval constraint `T_inf[secondary] > T_inf[source]` is enforced via a `-Inf` reject in the likelihood. The model is described in `METHODS.md`.

**Arguments**
- `d`: model data tuple as returned by [`build_data`](/api#TransmissionLinelist.build_data-Tuple{Any}), with fields `t0`, `onset_lo_day`, `onset_hi_day`, `exp_lo_day`, `exp_hi_day`, `source_idx`, `Zobs`, and `N`.
  
- `edges`: knot positions in days from `t0`, as returned by [`bin_edges_day`](/api#TransmissionLinelist.bin_edges_day-Tuple{Any}).
  

**Returns**

A `DynamicPPL.Model` ready to pass to `Turing.sample`. The sampled chain contains the population parameters `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`, `phi_inv_sqrt`, `σ_rw`, the derived `k` and `log_R`, the random-walk innovations `ε` and initial value `log_R_init`, and the per-case latent vectors `T_onset` and `T_inf`.

**Examples**

```julia
ll    = load_linelist()
d     = build_data(ll)
edges = bin_edges_day(d.t0)
m     = joint_model(d, edges)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/model.jl#L24" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.load_linelist' href='#TransmissionLinelist.load_linelist'><span class="jlbinding">TransmissionLinelist.load_linelist</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
load_linelist() -> DataFrames.DataFrame
load_linelist(path) -> DataFrames.DataFrame

```


Load and clean the Epuyén line list from a CSV file.

Reads the CSV at `path`, drops any duplicated rows with patient IDs ending in `_alt`, parses `exposure_lower`, `exposure_upper`, `onset_date`, `onset_lower`, and `onset_upper` as `Date`s (defaulting `onset_lower` and `onset_upper` to `onset_date` when absent), and sorts the rows by integer-valued `patient_id`. Used by [`build_data`](/api#TransmissionLinelist.build_data-Tuple{Any}) to produce the model input tuple.

**Arguments**
- `path`: path to a line-list CSV. Defaults to the bundled `data/linelist.csv` shipped with the package.
  

**Returns**

A `DataFrame` with one row per case and parsed date columns ready for [`build_data`](/api#TransmissionLinelist.build_data-Tuple{Any}).

**Examples**

```julia
using TransmissionLinelist
ll = load_linelist()
first(ll, 3)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/data.jl#L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.log_R_at-Tuple{Real, AbstractVector{<:Real}, Any}' href='#TransmissionLinelist.log_R_at-Tuple{Real, AbstractVector{<:Real}, Any}'><span class="jlbinding">TransmissionLinelist.log_R_at</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
log_R_at(
    t::Real,
    knots::AbstractVector{<:Real},
    log_R
) -> Any

```


Piecewise-linear interpolation of `log R(t)` between weekly knots.

Linearly interpolates `log_R` against `knots` at the time `t`, clamping to the endpoint values outside the knot range.

**Arguments**
- `t`: time (in days from `t0`) at which to evaluate log R.
  
- `knots`: knot positions in days, as returned by [`bin_edges_day`](/api#TransmissionLinelist.bin_edges_day-Tuple{Any}).
  
- `log_R`: vector of log R values at each knot.
  

**Returns**

The interpolated log R value at `t`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/data.jl#L150" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_data-Tuple{Any}' href='#TransmissionLinelist.plot_data-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.plot_data</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_data(ll)
```


Two-panel view of the raw line list: epicurve by ISO week of onset (left) and exposure windows against onset dates (right). Returns a `Makie.Figure`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/plots.jl#L13-L18" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_delta_sense_check-Tuple{Any, Any}' href='#TransmissionLinelist.plot_delta_sense_check-Tuple{Any, Any}'><span class="jlbinding">TransmissionLinelist.plot_delta_sense_check</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_delta_sense_check(chn, data)
```


Sense-check the per-pair posterior of δ against the fitted population `Normal(μ_δ, σ_δ)`. For each sourced pair, take the posterior of `δ_pair = T_inf[secondary] − T_onset[source]` and reduce to its median; then plot the histogram of those per-pair medians with the population density overlaid. Returns a `Makie.Figure`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/plots.jl#L384-L392" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_inc_sense_check-Tuple{Any, Any}' href='#TransmissionLinelist.plot_inc_sense_check-Tuple{Any, Any}'><span class="jlbinding">TransmissionLinelist.plot_inc_sense_check</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_inc_sense_check(chn, data; n_density_draws = 200)
```


Sense-check the per-case posterior of the incubation period against the fitted population `LogNormal(μ_inc, σ_inc)`. For each case, takes the posterior of `inc_i = T_onset[i] − T_inf[i]` and reduces to its median; plots the histogram of those per-case medians with the median PDF (and 95% pointwise ribbon) of the population LogNormal overlaid. Returns a `Makie.Figure`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/plots.jl#L618-L627" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_pair-Tuple{Any}' href='#TransmissionLinelist.plot_pair-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.plot_pair</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_pair(chn; thin = 2)
```


Corner plot of the population scalars `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`, `k` via PairPlots.jl. Returns a Makie `Figure` (requires a Makie backend such as CairoMakie loaded at the call site).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/plots.jl#L193-L199" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_predictive_distributions-Tuple{Any}' href='#TransmissionLinelist.plot_predictive_distributions-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.plot_predictive_distributions</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_predictive_distributions(chn; rng = Random.MersenneTwister(1))
```


Two-by-two panel of the implied population distributions under the posterior for incubation period, transmission timing δ, generation interval, and serial interval. Each panel shows draws from `p(y_new | data) = ∫ p(y_new | θ) p(θ | data) dθ`, i.e. what a new case or transmission pair would look like under the fitted parameters.

This is _not_ a posterior-predictive check against observed data; for that, see `plot_z_ppc`, `plot_delta_sense_check`, and `plot_inc_sense_check`.

Inc and δ panels overlay the parametric density (median PDF with a 95% pointwise ribbon across draws) and a histogram of one predictive sample per draw. GI and SI show the predictive-sample histogram only. Returns a `Makie.Figure`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/plots.jl#L266-L283" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_prior_predictives-Tuple{}' href='#TransmissionLinelist.plot_prior_predictives-Tuple{}'><span class="jlbinding">TransmissionLinelist.plot_prior_predictives</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_prior_predictives(; n = 5000, rng = Random.MersenneTwister(0))
```


Prior-predictive panel: histograms of Inc, δ, and GI/SI drawn from the package&#39;s independent priors on `μ_inc`, `σ_inc`, `μ_δ`, `σ_δ`. Returns a `Makie.Figure`.

Three histograms faceted by quantity is the kind of plot AoG was built for: one long-form data frame, `mapping(:value, layout = :panel)`, `visual(Hist)`. Each panel still has its own viewing window so long tails don&#39;t squash the bars; rather than per-facet axis limits, the input is pre-clipped to the window for each panel.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/plots.jl#L667-L679" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.plot_rt-Tuple{Any}' href='#TransmissionLinelist.plot_rt-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.plot_rt</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
plot_rt(chn; n_draws_plot = 100, ymax = 4.0)
```


Spaghetti plot of R(t) over the weekly knots. Each thinned posterior draw is a piecewise-linear trajectory through `(knot_date[b], exp(log_R[b]))`. Knot dates come from `BIN_EDGES` (data.jl). Returns a `Makie.Figure`.

Per-draw spaghetti is built as a long-form `DataFrame` and drawn via AlgebraOfGraphics with `group = :draw`, which is the idiomatic way to spell &quot;one line per draw&quot; once the data is tidy.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/plots.jl#L333-L343" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.prepare_model-Tuple{Any}' href='#TransmissionLinelist.prepare_model-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.prepare_model</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
prepare_model(
    ll
) -> Tuple{Any, NamedTuple{(:t0, :onset_lo_day, :onset_hi_day, :exp_lo_day, :exp_hi_day, :source_idx, :Zobs, :N), <:NTuple{8, Any}}, Vector{Float64}}

```


Build the joint model from a line-list `ll`.

Wraps [`build_data`](/api#TransmissionLinelist.build_data-Tuple{Any}), [`bin_edges_day`](/api#TransmissionLinelist.bin_edges_day-Tuple{Any}), and [`joint_model`](/api#TransmissionLinelist.joint_model-Tuple{Any,%20Any}) into a single call so the analysis walkthrough and CLI share the same model construction code path.

**Arguments**
- `ll`: a line-list `DataFrame` as returned by [`load_linelist`](/api#TransmissionLinelist.load_linelist).
  

**Returns**

A 3-tuple `(model, d, edges)`: the Turing model, the augmented data named tuple from [`build_data`](/api#TransmissionLinelist.build_data-Tuple{Any}), and the weekly knot edges from [`bin_edges_day`](/api#TransmissionLinelist.bin_edges_day-Tuple{Any}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/data.jl#L124" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.sample_fit-Tuple{Any}' href='#TransmissionLinelist.sample_fit-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.sample_fit</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sample_fit(model; samples=1000, chains=4, target_accept=0.95,
           seed=20260508, progress=false)
```


Run NUTS on `model` using the package&#39;s default Enzyme AD backend and `InitFromPrior()` chain initialisation. Returns the FlexiChain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/main.jl#L4-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.save_posterior-Tuple{Any, Any}' href='#TransmissionLinelist.save_posterior-Tuple{Any, Any}'><span class="jlbinding">TransmissionLinelist.save_posterior</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
save_posterior(post, path) -> Any

```


Write the posterior summary `post` (as returned by [`summarise`](/api#TransmissionLinelist.summarise-Tuple{Any})) to a CSV at `path`, one column per scalar parameter plus one column per `log_R` knot.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/postprocess.jl#L106" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.summarise-Tuple{Any}' href='#TransmissionLinelist.summarise-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.summarise</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
summarise(
    chn
) -> NamedTuple{(:μ_inc, :σ_inc, :μ_δ, :σ_δ, :k, :log_R_chain, :mean_gi_si, :sd_gi_si, :p_pre), <:Tuple{Any, Any, Any, Any, Any, Any, Any, Any, Dict{Float64, Vector{Float64}}}}

```


Build the named-tuple of posterior draws consumed by [`save_posterior`](/api#TransmissionLinelist.save_posterior-Tuple{Any,%20Any}) and print the headline summary table via [`summary_table`](/api#TransmissionLinelist.summary_table-Tuple{Any}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/postprocess.jl#L72" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.summary_table-Tuple{Any}' href='#TransmissionLinelist.summary_table-Tuple{Any}'><span class="jlbinding">TransmissionLinelist.summary_table</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
summary_table(chn)
```


Posterior summary `DataFrame` for the headline quantities: incubation mean, 95th and 99th percentiles, transmission timing μ_δ / σ_δ, GI / SI mean and SD, and Negative-Binomial dispersion k. Columns: `quantity`, `median`, `lower_95`, `upper_95`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/plots.jl#L123-L130" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.vector_chain-Tuple{Any, Symbol}' href='#TransmissionLinelist.vector_chain-Tuple{Any, Symbol}'><span class="jlbinding">TransmissionLinelist.vector_chain</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
vector_chain(chn, name::Symbol) -> Any

```


Return a vector of pooled posterior samples for each entry of a vector-valued parameter (e.g. `:T_inf`, `:log_R`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/postprocess.jl#L49" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TransmissionLinelist.z_ppc_summary-Tuple{Any, Any}' href='#TransmissionLinelist.z_ppc_summary-Tuple{Any, Any}'><span class="jlbinding">TransmissionLinelist.z_ppc_summary</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
z_ppc_summary(chn, d; rng = Random.MersenneTwister(1),
              edges = bin_edges_day(d.t0))
```


Companion to `plot_z_ppc` returning a `DataFrame` of numeric posterior-predictive summaries for three discrete test statistics — `sum(Z)`, `max(Z)`, and `count(Z = 0)`. Replicates `Z_rep` jointly in `(T_inf, log_R, k)` to match `plot_z_ppc`. Columns: `statistic`, `observed`, `rep_median`, `rep_lower_95`, `rep_upper_95`, `p_ppp`, where `p_ppp = 2 · min(P(T_rep ≥ T_obs), P(T_rep ≤ T_obs))` is the two-sided Bayesian posterior-predictive p-value.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/epiforecasts/andv-linelist-analysis/blob/c482faa10bcfe855bf238d25ef40e503f4f71a6a/src/plots.jl#L494-L505" target="_blank" rel="noreferrer">source</a></Badge>

</details>

