# Limitations

Known caveats of the joint model fit to the Epuyén 2018–19 Andes hantavirus outbreak.
See [MODEL.md](https://github.com/sbfnk/hantavirus/blob/main/MODEL.md) for the model description.

## Exposure encoding pins transmission-timing variability

Most of what the model can say about transmission timing is limited by how the line list was recorded.
31 of 33 sourced pairs have a single-day exposure window, and that day is almost always the source's symptom onset.
The fitted transmission-timing SD of about 0.6 d mostly reflects within-day uncertainty in `T_inf` rather than biological spread.
The two cannot be disentangled from these data.
Multi-day pre-symptomatic transmission is therefore rare in this outbreak (P(δ < −1 d) ≈ 3%, P(δ < −2 d) essentially zero).
The split into "any pre-symptomatic" vs "post-symptomatic" would be dominated by this within-day floor and is not reported.

## Late R(t) bins are prior-driven

There are very few cases after early January 2019, and the random walk on `log R(t)` reverts to its prior in those bins.
The wide credible intervals on the right of the R(t) figure show this.

## Offspring dispersion `k` has prior dependence

34 cases is thin for identifying a Negative-Binomial dispersion.
The prior on `k`, centred at 0.3, has visible influence on the posterior centre.

## Right-truncation of long incubation periods

The Martínez paper does not document a surveillance cut-off date.
The last observed onset is 2019-02-06.
If surveillance effectively ended shortly after, a case infected close to the cut-off with a long incubation period could have been missed because the onset would land after the cut-off.
This biases the upper tail of the incubation distribution slightly downward, affecting mainly the ~3–5 late-infected cases (those infected from late January 2019 onward).
Patient 1 and the bulk of the line list were infected weeks before the last onset, so long incubation periods would have had time to materialise.
The reported 99th percentile of about 45 d is therefore a mild lower bound.

## Offspring count Z is restricted to high-certainty transmissions

The paper notes that *"only events of person-to-person transmission with a high certainty of infection at the time of the event were included in analyses and are reported as Z values."*
Weakly-attributed transmissions are dropped from Z, so the observed Z is a lower bound on true offspring count per case.
This biases R(t) somewhat downward.
The effect on Negative-Binomial dispersion `k` depends on whether dropped events were concentrated on high-Z cases (would push `k` higher, toward Poisson) or spread evenly (would push `k` lower).

## Citing reporting practices

The reporting structure of this document follows:

> Charniga K, et al. *Best practices for estimating and reporting epidemiological delay distributions of infectious diseases.* 2024. [arXiv:2405.08841](https://arxiv.org/abs/2405.08841)
