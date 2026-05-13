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
