# Limitations

Known caveats of the joint model fit to the Epuyén 2018–19 Andes hantavirus outbreak.
See [MODEL.md](https://github.com/sbfnk/hantavirus/blob/main/MODEL.md) for the model description.

## Right-truncation of long incubation periods

The Martínez paper does not document a surveillance cut-off date.
The last observed onset is 2019-02-06.
If surveillance effectively ended shortly after, a case infected close to the cut-off with a long incubation period could have been missed because the onset would land after the cut-off.
This biases the upper tail of the incubation distribution slightly downward, affecting mainly the ~3–5 late-infected cases (those infected from late January 2019 onward).
Patient 1 and the bulk of the line list were infected weeks before the last onset, so long incubation periods would have had time to materialise.
The reported 99th percentile of about 45 d is therefore a mild lower bound.

## Offspring dispersion `k` has prior dependence

34 cases is thin for identifying a Negative-Binomial dispersion.
The prior on `k`, centred at 0.3, has visible influence on the posterior centre.

## Offspring count Z is restricted to high-certainty transmissions

The paper notes that *"only events of person-to-person transmission with a high certainty of infection at the time of the event were included in analyses and are reported as Z values."*
Weakly-attributed transmissions are dropped from Z, so the observed Z is a lower bound on true offspring count per case.
This biases R(t) somewhat downward.
The effect on Negative-Binomial dispersion `k` depends on whether dropped events were concentrated on high-Z cases (would push `k` higher, toward Poisson) or spread evenly (would push `k` lower).

## Late R(t) bins are prior-driven

There are very few cases after early January 2019, and the random walk on `log R(t)` reverts to its prior in those bins.
The wide credible intervals on the right of the R(t) figure show this.

## Real-time fitting caveats

The real-time corrections handle three specific biases — long-incubation cases, late transmissions, and incomplete clusters.
Not corrected:

- geographic / severity / surveillance reporting biases,
- the onset-to-report delay (only chain completion is modelled),
- general under-ascertainment,
- incomplete source attribution,
- pre-symptomatic transmission with an unobserved source,
- ongoing zoonosis.

**Pre-symptomatic transmission with an unobserved source.**
When `δ < −Inc[src]`, a source's onset can be later than its secondary's onset.
At an `obs_time` cut-off the secondary can be in the line list while the source isn't; `filter_realtime` then drops the source attribution and the secondary looks like an apparent index.
Probably small for ANDV (δ averages near zero with σ_δ ≈ 1) but a real selection effect the current implementation doesn't correct for.

**Ongoing zoonosis.**
The model treats index (zoonotic) cases as a small starter set; cluster-completeness only thins observed sources, it doesn't add back population members whose Inc hasn't completed yet.
The current implementation is fine for an outbreak with a few initial zoonotic cases and no ongoing zoonosis (the Epuyén pattern); it would under-count cases if zoonosis continued throughout the outbreak.
