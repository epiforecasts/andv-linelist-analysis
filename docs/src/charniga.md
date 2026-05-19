# Reporting best practices (Charniga et al.)

This page maps the joint-model reporting in this package against the
checklist in:

> Charniga K, Park SW, Akhmetzhanov AR, Cori A, Dushoff J, Funk S,
> Gostic KM, Linton NM, Lison A, Overton CE, Pulliam JRC, Ward T,
> Cauchemez S, Abbott S. *Best practices for estimating and reporting
> epidemiological delay distributions of infectious diseases.* PLOS
> Computational Biology, 2024.
> [doi:10.1371/journal.pcbi.1012520][doi]

[doi]: https://doi.org/10.1371/journal.pcbi.1012520

The Charniga et al. checklist (their Table 2) is split into items
covering estimation and items covering reporting. The rows below
paraphrase each item. Wording is condensed; consult the paper for the
authoritative version.

For items marked "Not yet" a tracking issue is referenced where one
exists. Items marked with `(?)` are uncertain and need a closer look at
the published table.

## Estimation

| # | Recommendation (paraphrased) | Addressed? | Where in this package |
|---|---|---|---|
| E1 | Adjust for interval censoring of event times | ✓ | [Model — Latent infection and onset times](model.md). Each case has continuous latents `T_inf[i]`, `T_onset[i]` sampled uniformly inside the recorded exposure and onset windows; the likelihood is evaluated at the latents, not at window midpoints. Implemented in `joint_model` in `src/model.jl`. |
| E2 | Adjust for right truncation when fitting on data collected mid-outbreak | ✓ | [Real-time monitoring](realtime.md) and `truncation_model` in `src/model.jl`. Per-case `-logcdf(inc_dist, obs_time - T_inf[i])` plus an offspring-completeness factor through `ConvolvedDelays`. |
| E3 | Adjust for dynamical / epidemic-phase bias in delay estimation | ✓ (by construction) | Delays are fit jointly with `R(t)` rather than from forward-time pair observations, so the epidemic-phase reweighting is handled inside the joint likelihood rather than by a post-hoc correction. Discussed in [Model](model.md). |
| E4 | Fit multiple candidate probability distributions and compare them | Partial | Only LogNormal is fitted for the incubation period and Normal for transmission timing δ. No formal model comparison (WAIC, LOO) across candidate families is run. Tracked as **proposed issue: "Fit alternative distributional families (Gamma, Weibull) for the incubation period and compare via PSIS-LOO"**. |
| E5 | Convert between parameterisations carefully (mean/SD vs shape/scale etc.) | ✓ | Conversions are done through `Distributions.jl` constructors (`LogNormal(μ, σ)`, `Normal(μ, σ)`) and `quantile(LogNormal(...), p)` rather than hand-rolled formulae. See `src/plots.jl` (`summary_table`) and `MODEL.md`. |
| E6 | Visualise the fitted distribution alongside the observed data and observation process | ✓ | [Analysis walkthrough](analysis.md) posterior-predictive sections render the fitted incubation density, transmission-timing density, and offspring distribution against observed pairs and Z counts. Built from `plot_inc_posterior_predictive`, `plot_delta_posterior_predictive`, and `plot_rep_posterior_predictive` in `src/plots.jl`. |
| E7 | Stratify estimates by relevant subgroups when sample size permits | Not yet | 34 cases is too thin to stratify reliably (no age/sex/severity stratification is run). The line list itself does not record those covariates. Tracked as **proposed issue: "Document why stratified estimates are not produced for the Epuyén line list"**. |
| E8 | Check Bayesian model diagnostics (R̂, ESS, divergences, trace plots) | ✓ | `diagnostics_table` in `src/plots.jl` and the *Diagnostics* section of [Analysis walkthrough](analysis.md) report maximum R̂, minimum bulk ESS, divergence count, and wall-clock sampling time. Trace plots are not currently rendered in the docs. Tracked as **proposed issue: "Add trace/rank plots for key parameters to the diagnostics section"**. |

## Reporting

| # | Recommendation (paraphrased) | Addressed? | Where in this package |
|---|---|---|---|
| R1 | Report a measure of variability (SD / dispersion), not only a central tendency | ✓ | `summary_table` in `src/plots.jl` reports posterior medians with 95% CrIs for `μ`, `σ`, the SD of the fitted distributions, and offspring dispersion `k`. The README "Headline results" table mirrors this. |
| R2 | Report quantiles of the fitted distribution (e.g. 95th, 99th percentile) | ✓ | `summary_table` reports the 95th and 99th percentiles of the fitted incubation period with 95% CrIs. See *Key outputs* in [Analysis walkthrough](analysis.md). |
| R3 | State the probability density function / family explicitly | ✓ | [Model](model.md) writes out `Inc ∼ LogNormal(μ_inc, σ_inc)`, `δ ∼ Normal(μ_δ, σ_δ)`, and the offspring `NegativeBinomial(R, k)` parameterisation. |
| R4 | Report credible or confidence intervals (typically 90% or 95%) | ✓ | 95% credible intervals throughout. The interval level is stated in the README and [Analysis walkthrough](analysis.md). |
| R5 | Include contextual information (sample size, epidemic curve, control measures, location, time) | Partial | Sample size (34 cases), location (Epuyén, Argentina), and outbreak window (2018–19) are reported in the README, [index](index.md), and [Analysis walkthrough](analysis.md). The epidemic curve is rendered by `plot_data`. Public-health control measures (cordon, contact tracing) are referenced via the Martínez et al. citation but not summarised in the docs. Tracked as **proposed issue: "Add a short context paragraph on Epuyén outbreak control measures with citations"**. |
| R6 | Report demographic and exposure-route summaries of the contributing cases | Not yet | The line list does not carry age/sex/exposure-route covariates beyond `source_case` attribution, so no demographic breakdown is shown. Tracked as **proposed issue: "Document demographic/exposure-route information available in the source publication and add a static summary table"**. |
| R7 | Share code and (anonymised) data in an open repository | ✓ | Code is on GitHub at `epiforecasts/andv-linelist-analysis`; the hand-encoded line list is shipped in `data/linelist.csv` and loaded by `load_linelist`. License is in `LICENSE`. |
| R8 | Provide the line-list-format data with the boundaries required for re-estimation (exposure / onset windows, attributed source) | ✓ | `data/linelist.csv` contains `exposure_lower`, `exposure_upper`, `onset_date` (with the `_alt` sensitivity rows), and `source_case` / `Z` columns. The parser `load_linelist` and the data-construction step `build_data` in `src/data.jl` are documented in the [API Reference](api.md). |
| R9 | Document known limitations / sources of residual bias | ✓ | A dedicated [Limitations page](limitations.md) covers right-truncation of long incubation periods, prior dependence of offspring dispersion `k`, the high-certainty filter on `Z`, prior-driven late `R(t)` bins, and real-time-fitting caveats (geography, severity, reporting delay, under-ascertainment, pre-symptomatic transmission with an unobserved source, ongoing zoonosis). |

## Summary

Of the items extracted from the Charniga et al. checklist, the joint
model and its documentation cover most directly. The main gaps are
multi-distribution model comparison (E4), stratification by subgroup
(E7), trace/rank plots (E8), an explicit context paragraph on control
measures (R5), and demographic/exposure-route summaries (R6). Each of
these has a proposed tracking issue listed above.
