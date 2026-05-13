# Model

The joint Turing model for the TransmissionLinelist.jl analysis of the Epuyén 2018–19 Andes hantavirus outbreak.
See the [README](https://github.com/sbfnk/hantavirus/blob/main/README.md) for headline results and [LIMITATIONS.md](https://github.com/sbfnk/hantavirus/blob/main/LIMITATIONS.md) for known caveats.

Each case has two continuous latents.
`T_onset[i]` is uniform over the recorded onset window, which is one day wide if only a single onset date was recorded.
`T_inf[i]` is uniform over the exposure window for sourced cases, or over an 80-day pre-onset window for the zoonotic index.

| Quantity | Distribution | Priors |
|---|---|---|
| Incubation period (`T_onset − T_inf`) | LogNormal | log-mean ~ Normal(3.0, 0.5), log-SD ~ half-Normal(0, 0.5) |
| Transmission timing relative to source onset (`T_inf(sec) − T_onset(src)`) | Normal | mean ~ Normal(0, 5), SD ~ half-Normal(0, 1) |
| Offspring count `Z` per case | Negative-Binomial with mean `R(t)` and dispersion `k` | `k` ~ half-Normal(0.3, 0.5) |
| `log R(t)` at weekly knots, linearly interpolated between them | Random walk | first knot ~ Normal(log 1.5, 1); innovation SD ~ half-Normal(0, 0.5) |

A per-pair constraint enforces `T_inf(secondary) > T_inf(source)` so the generation interval is positive.
Generation interval is transmission timing plus the source's incubation period.
Serial interval is transmission timing plus the secondary's incubation period.
Both are derived in post-processing.

Inference uses NUTS, 4 chains, 1000 post-warmup samples each, `target_accept = 0.95`.
Default seed: 20260508.
