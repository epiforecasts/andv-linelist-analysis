# Model

The model estimates four quantities jointly from the Epuyén line list: the incubation period of Andes hantavirus, the timing of each documented onward transmission relative to its source case's symptom onset, a weekly time-varying reproduction number R(t), and the over-dispersion of offspring counts per case.

Because exposure and symptom-onset dates in the line list are recorded as windows rather than exact dates, each case is given two latent variables: a true infection time `T_inf` and a true symptom-onset time `T_onset`.
The model jointly explores these latent times together with the parameters governing the four quantities above.

Generation and serial intervals are not given separate priors.
They are computed from the fitted incubation period and transmission timing distributions in post-processing: the generation interval is the transmission timing plus the source's incubation period, and the serial interval is the transmission timing plus the secondary's incubation period.

See the [README](https://github.com/epiforecasts/andv-linelist-analysis/blob/main/README.md) for headline results and [LIMITATIONS.md](https://github.com/epiforecasts/andv-linelist-analysis/blob/main/LIMITATIONS.md) for known caveats.

## Reproduction number

The `R(t)` reported here is the case reproduction number indexed by source symptom onset: the expected number of secondary infections produced by a case whose onset is at time `t`.
Onset is the natural choice for these data: fitted transmission timing is tightly clustered around source onset (`μ_δ ≈ 0`, `σ_δ ≈ 0.6 d`), so a case's offspring are infected within roughly a day of the case becoming symptomatic.
Because transmission is so concentrated at onset, this mostly coincides with the instantaneous reproduction number indexed by infection date.

## Latent variables and priors

Each case has two continuous latent variables.
`T_onset[i]` has a uniform prior over the recorded onset window, which is one day wide if only a single onset date was recorded.
`T_inf[i]` has a uniform prior over the exposure window for sourced cases, or over an 80-day pre-onset window for the zoonotic index.

| Quantity | Distribution | Priors |
|---|---|---|
| Incubation period (`T_onset − T_inf`) | LogNormal | log-mean ~ Normal(3.0, 0.5), log-SD ~ half-Normal(0, 0.5) |
| Transmission timing relative to source onset (`T_inf(sec) − T_onset(src)`) | Normal | mean ~ Normal(0, 5), SD ~ half-Normal(0, 1) |
| Offspring count `Z` per case | Negative-Binomial with mean `R(t)` and dispersion `k` | `k` ~ half-Normal(0.3, 0.5) |
| `log R(t)` at weekly knots, linearly interpolated between them | Random walk | first knot ~ Normal(log 1.5, 1); innovation SD ~ half-Normal(0, 0.5) |

A per-pair constraint enforces `T_inf(secondary) > T_inf(source)` so the generation interval is positive.

Inference uses NUTS, 4 chains, 1000 post-warmup samples each, `target_accept = 0.95`.
Default seed: 20260508.
