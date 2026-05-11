# Andes virus — joint estimation of incubation, transmission timing, and R(t)

A Julia + Turing model for the Epuyén 2018–19 Andes hantavirus outbreak
([Martínez et al. 2020, NEJM](https://doi.org/10.1056/NEJMoa2009040)).
From the line list given in the paper it jointly estimates the incubation
period and the transmission timing of each secondary infection relative to
its source's symptom onset, plus a time-varying reproduction number with
offspring dispersion. Double interval censoring of exposure and onset is
handled by Bayesian data augmentation over continuous latent infection and
onset times.

The generation interval (transmission timing plus the source's incubation
period) and the serial interval (transmission timing plus the secondary's
incubation period) are derived in post-processing from the fitted
distributions. A per-pair constraint that the secondary's infection time is
later than the source's keeps the generation interval positive at the
latent level.

## Model

For each case `i` the model has continuous latents:

- `T_onset[i]` ~ Uniform over the recorded onset window (defaults to a
  one-day window when only a single onset date was recorded).
- `T_inf[i]` ~ Uniform over the exposure window (sourced cases) or a
  wide pre-onset window of 80 days (the zoonotic index case).

The estimated quantities and their distributional choices:

| Quantity | Distribution | Parameters |
|---|---|---|
| Incubation period `Inc = T_onset − T_inf` | LogNormal(μ_inc, σ_inc) | `μ_inc ~ Normal(3.0, 0.5)`, `σ_inc ~ half-Normal(0, 0.5)` |
| Transmission timing `δ = T_inf(sec) − T_onset(src)` | Normal(μ_δ, σ_δ) | `μ_δ ~ Normal(0, 5)`, `σ_δ ~ half-Normal(0, 1)` |
| Offspring count `Z` (per case) | Negative-Binomial(`k`, `k/(k + R(t))`) | `k ~ half-Normal(0.3, 0.5)`, truncated at 0 |
| Log reproduction number `log R(t)` | random walk with innovation SD `σ_rw` over weekly bins | `log R[1] ~ Normal(log 1.5, 1)`, `σ_rw ~ half-Normal(0, 0.5)` |

Generation interval (`GI = δ + Inc(source)`) and serial interval
(`SI = δ + Inc(secondary)`) are derived from posterior draws.

Inference is by NUTS (4 chains × 1000 samples after warmup,
`target_accept = 0.95`). Reproducibility seed is set in `scripts/run.jl`.

## Headline results (Epuyén line list)

### Incubation period

Distribution: **LogNormal**.

| Quantity | Posterior median (95% CrI) |
|---|---|
| Mean | 22.6 d (20.2 – 25.5) |
| 95th percentile | 36.2 d (31.3 – 44.5) |
| 99th percentile | 45.0 d (37.6 – 58.7) |

### Transmission timing relative to source onset

Distribution: **Normal**. Negative values mean the secondary was infected before the source became symptomatic.

| Quantity | Posterior median (95% CrI) |
|---|---|
| Mean | 0.2 d (−0.2 – 0.5) |
| SD | 0.6 d (0.5 – 0.8) |
| P(transmission is pre-symptomatic) | 0.4 (0.2 – 0.6) |
| P(pre-symptomatic by more than 1 day) | 0.03 (0.00 – 0.12) |
| P(pre-symptomatic by more than 2 days) | 0.00 (0.00 – 0.01) |

### Generation interval / serial interval

Derived from incubation and transmission timing.

| Quantity | Posterior median (95% CrI) |
|---|---|
| Mean | 22.7 d (20.4 – 25.6) |
| SD | 7.4 d (5.7 – 10.6) |

### Offspring distribution

Distribution: **Negative-Binomial** with mean `R(t)` and dispersion `k`.

| Quantity | Posterior median (95% CrI) |
|---|---|
| Dispersion `k` | 0.4 (0.1 – 0.9) |

### Time-varying reproduction number R(t)

Weekly bins. Late-outbreak bins fall back toward the random-walk prior (few cases).

| Week | R median (95% CrI) |
|---|---|
| ≤ 2018-11-12 | 1.64 (0.65 – 4.87) |
| 2018-11-12 – 2018-11-19 | 1.50 (0.48 – 5.97) |
| 2018-11-19 – 2018-11-26 | 1.43 (0.59 – 4.01) |
| 2018-11-26 – 2018-12-03 | 0.96 (0.27 – 3.39) |
| 2018-12-03 – 2018-12-10 | 0.63 (0.15 – 2.30) |
| 2018-12-10 – 2018-12-17 | 0.40 (0.12 – 1.49) |
| 2018-12-17 – 2018-12-24 | 0.31 (0.04 – 1.63) |
| 2018-12-24 – 2018-12-31 | 0.24 (0.02 – 1.58) |
| 2018-12-31 – 2019-01-07 | 0.20 (0.01 – 1.65) |
| 2019-01-07 – 2019-01-14 | 0.20 (0.01 – 1.75) |
| 2019-01-14 – 2019-01-21 | 0.19 (0.01 – 2.14) |
| 2019-01-21 – 2019-01-28 | 0.20 (0.00 – 2.40) |
| 2019-01-28 – 2019-02-04 | 0.20 (0.00 – 3.36) |
| > 2019-02-04 | 0.20 (0.00 – 3.76) |

## Limitations

- **σ_δ is dominated by within-day exposure encoding rather than biological variation.**
  In the Martínez line list, 31 of 33 sourced pairs have an exposure window
  equal to a single day, most often the source's symptom onset day (the
  documented contact event). The fitted σ_δ ≈ 0.6 d therefore reflects
  within-day uncertainty in `T_inf` rather than the true biological spread of
  transmission timing. The reported pre-symptomatic transmission fraction is
  conditional on the assumption that the recorded contact event is when
  transmission occurred; if actual contacts spanned days, pre-symptomatic
  transmission may be under-recorded.
- **Late-outbreak R(t) bins are prior-driven.** Few cases occur after early
  January 2019, so the random walk falls back toward its prior (visible as
  wide credible intervals).
- **Offspring dispersion `k` has some prior dependence.** 34 cases is thin for
  identifying a Negative-Binomial dispersion; the prior centred at 0.3
  visibly nudges the posterior centre.

## Repository layout

```
src/
  data.jl          — line list loading and bin definitions
  model.jl         — the joint Turing model (incubation, transmission timing, R(t))
  postprocess.jl   — diagnostics, summaries, CSV output
scripts/
  run.jl           — entry point
data/
  linelist.csv     — Epuyén outbreak line list (Martínez Table S2)
Project.toml       — Julia environment
LICENSE            — MIT
```

## Data

The Epuyén line list (`data/linelist.csv`) is hand-encoded from Table S2 of
the supplementary appendix of Martínez et al. 2020. Columns: patient ID, age,
sex, residence, exposure place, exposure window (lower / upper), onset date,
attributed source (or `index` for the zoonotic case), relationship to source,
transmission wave, observed offspring count `Z`, and free-text notes.

## Running

```
julia --project=. -t auto scripts/run.jl
```

NUTS, 4 chains × 1000 samples. Takes a few minutes on a laptop. Posterior
saved to `output/posterior.csv`.

## Citing

If you use this code or the Epuyén line list encoding, please cite:

> Martínez VP, Di Paola N, Alonso DO, et al. *"Super-spreaders" and
> person-to-person transmission of Andes virus in Argentina.* N Engl J Med
> 2020;383:2230–41. [doi:10.1056/NEJMoa2009040](https://doi.org/10.1056/NEJMoa2009040)

We follow the reporting recommendations of:

> Charniga K, et al. *Best practices for estimating and reporting
> epidemiological delay distributions of infectious diseases.* 2024.
> [arXiv:2405.08841](https://arxiv.org/abs/2405.08841)

## License

MIT (see [LICENSE](LICENSE)).
