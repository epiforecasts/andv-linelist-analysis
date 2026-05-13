# Andes virus — Epuyén 2018–19 line-list analysis

A Julia + Turing joint model fitted to the line list of the
[Epuyén 2018–19 ANDV outbreak](https://doi.org/10.1056/NEJMoa2009040).
Estimates the incubation period, transmission timing of each secondary
relative to its source's symptom onset, a weekly time-varying reproduction
number, and offspring dispersion jointly from interval-censored exposure
and onset windows.

## Pages

- [Methods](methods.md) — priors, data augmentation, GI / SI derivation, limitations.
- [API Reference](api.md) — exported functions.

Headline posterior numbers, repository layout, and run instructions live in
the [README on GitHub](https://github.com/epiforecasts/andv-linelist-analysis).
