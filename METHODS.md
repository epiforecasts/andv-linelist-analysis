# Methods and limitations

Detailed model description, prior choices, and known limitations for the
TransmissionLinelist.jl analysis of the Epuyén 2018–19 Andes hantavirus outbreak. See
the [README](README.md) for headline results and how to run the analysis.

## Model

Each case has two continuous latents. `T_onset[i]` is uniform over the
recorded onset window, which is one day wide if only a single onset date
was recorded. `T_inf[i]` is uniform over the exposure window for sourced
cases, or over an 80-day pre-onset window for the zoonotic index.

| Quantity | Distribution | Priors |
|---|---|---|
| Incubation period (`T_onset − T_inf`) | LogNormal | log-mean ~ Normal(3.0, 0.5), log-SD ~ half-Normal(0, 0.5) |
| Transmission timing relative to source onset (`T_inf(sec) − T_onset(src)`) | Normal | mean ~ Normal(0, 5), SD ~ half-Normal(0, 1) |
| Offspring count `Z` per case | Negative-Binomial with mean `R(t)` and dispersion `k` | `1/√k` ~ half-Normal(0, 1) (Stan's reciprocal-sqrt parameterisation; diffuse and symmetric on the overdispersion SD-multiplier scale) |
| `log R(t)` over weekly bins | Non-centred random walk | first bin ~ Normal(log 1.5, 1); innovation SD ~ half-Normal(0, 0.2) (≈ 5%/day typical, ≈ 15%/day at the prior 95th percentile) |

A per-pair constraint enforces `T_inf(secondary) > T_inf(source)` so that
the generation interval is positive. Generation interval = transmission
timing + source's incubation period; serial interval = transmission timing
+ secondary's incubation period. Both are computed in post-processing.

Inference uses NUTS, 4 chains, 1000 post-warmup samples each, `target_accept = 0.95`. Default seed: 20260508.

## Limitations

### Exposure encoding pins transmission-timing variability

Most of what the model can say about transmission timing is limited by how
the line list was recorded. 31 of 33 sourced pairs have a single-day
exposure window, and that day is almost always the source's symptom onset.
The fitted transmission-timing SD of about 0.6 d mostly reflects within-day
uncertainty in `T_inf` rather than biological spread; we cannot disentangle
the two from these data. Multi-day pre-symptomatic transmission is therefore
robustly rare in this outbreak (P(δ < −1 d) ≈ 3%, P(δ < −2 d) essentially
zero), but the split into "any pre-symptomatic" vs "post-symptomatic" would
be dominated by this within-day floor and is not reported.

### Late R(t) bins are prior-driven

There are very few cases after early January 2019, and the random walk on
`log R(t)` reverts to its prior in those bins. The wide credible intervals
on the right of the R(t) figure show this.

### Offspring dispersion `k` is data-thin

34 cases is thin for identifying a Negative-Binomial dispersion. The
posterior on `k` is wide; the reciprocal-sqrt prior on `1/√k` is diffuse
and symmetric on the overdispersion scale so it does not pull the
posterior toward any particular value, but the credible interval is
necessarily wide given the sample size.

### Right-truncation of long incubation periods

The Martínez paper does not document a surveillance cut-off date. The last
observed onset is 2019-02-06. If surveillance effectively ended shortly
after, a case infected close to the cut-off with a long incubation period
could have been missed (onset would land after the cut-off). This biases
the upper tail of our incubation distribution slightly downward, affecting
mainly the ~3–5 late-infected cases (those infected from late January 2019
onward). Patient 1 and the bulk of the line list were infected weeks before
the last onset, so long incubation periods would have had time to
materialise. The reported 99th percentile of about 45 d is therefore a
mild lower bound.

### Offspring count Z is restricted to high-certainty transmissions

The paper notes that *"only events of person-to-person transmission with a
high certainty of infection at the time of the event were included in
analyses and are reported as Z values."* Weakly-attributed transmissions
are dropped from Z, so the observed Z is a lower bound on true offspring
count per case. This biases R(t) somewhat downward; the effect on Negative-
Binomial dispersion `k` depends on whether dropped events were concentrated
on high-Z cases (would push `k` higher, i.e., toward Poisson) or spread
evenly (would push `k` lower).

## Citing reporting practices

The reporting structure of this document follows:

> Charniga K, et al. *Best practices for estimating and reporting
> epidemiological delay distributions of infectious diseases.* 2024.
> [arXiv:2405.08841](https://arxiv.org/abs/2405.08841)
