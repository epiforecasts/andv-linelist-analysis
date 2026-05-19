# Andes virus — joint estimation of incubation, transmission timing, and R(t)

[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://epiforecasts.io/andv-linelist-analysis/dev)

Full documentation — model description, analysis walkthrough, API reference — at
<https://epiforecasts.io/andv-linelist-analysis/dev>.

A Julia + Turing model fitted to the Epuyén 2018–19 Andes hantavirus outbreak
([Martínez et al. 2020, NEJM](https://doi.org/10.1056/NEJMoa2009040)).

The model estimates four things from the line list in the paper: the
incubation period, the transmission timing of each secondary infection
relative to its source's symptom onset, a weekly time-varying reproduction
number, and offspring dispersion. Exposure and onset dates are
interval-censored. The model handles that by giving each case a continuous
latent infection time and a continuous latent onset time, each sampled
within its recorded window. Generation interval and serial interval are
derived from the fitted distributions in post-processing.

## Results

Rendered walkthrough with all tables and figures regenerated from the current model:
<https://epiforecasts.io/andv-linelist-analysis/dev/analysis>.

Raw artefacts from the most recent main build (`output/posterior.csv` and all figures):
<https://github.com/epiforecasts/andv-linelist-analysis/releases/tag/main-latest>.

## Methods and limitations

Model description and priors are in [MODEL.md](MODEL.md).
Known caveats are in [LIMITATIONS.md](LIMITATIONS.md).

## Repository layout

```
src/
  TransmissionLinelist.jl    — module entry point and imports
  data.jl          — line list loading and bin definitions
  model.jl         — the joint Turing model (incubation, transmission timing, R(t))
  plots.jl         — figure construction
  postprocess.jl   — diagnostics, summaries, CSV output
  main.jl          — CLI entry point (argument parsing)
data/
  linelist.csv     — Epuyén outbreak line list (Martínez Table S2)
docs/              — Documenter site with the analysis walkthrough
Project.toml       — Julia package manifest
Manifest.toml      — locked dependency versions
LICENSE            — MIT
```

Posterior and figures are regenerated locally by `analyse()` and published
to the [`main-latest`](https://github.com/epiforecasts/andv-linelist-analysis/releases/tag/main-latest)
release on every push to `main`; neither is committed.

## Data

`data/linelist.csv` is hand-encoded from Table S2 of the supplementary
appendix of Martínez et al. 2020. Columns: patient ID, age, sex, residence,
exposure place, exposure window (lower / upper), onset date, attributed
source (or `index` for the zoonotic case), relationship to source,
transmission wave, observed offspring count `Z`, and free-text notes.

## Running

```
julia --project=. -t auto -m TransmissionLinelist
```

A few minutes on a laptop. Posterior saved to `output/posterior.csv` and
figures to `figures/`.

Options:

```
-d, --data      path to linelist CSV   (default: data/linelist.csv)
-o, --output    output directory        (default: output/)
-f, --figures   figures directory       (default: figures/)
-n, --samples   NUTS samples per chain  (default: 1000)
-c, --chains    number of chains        (default: 4)
-s, --seed      random seed             (default: 20260508)
```

Example:

```
julia --project=. -t auto -m TransmissionLinelist -- -n 500 -c 2 -o results/
```

### From the REPL

```julia
julia> using TransmissionLinelist
julia> analyse()                                           # all defaults
julia> analyse(chains=2, samples=500, output="results/")  # with options
```

## Citing

If you use this code or the line list encoding, please cite:

> Martínez VP, Di Paola N, Alonso DO, et al. *"Super-spreaders" and
> person-to-person transmission of Andes virus in Argentina.* N Engl J Med
> 2020;383:2230–41. [doi:10.1056/NEJMoa2009040](https://doi.org/10.1056/NEJMoa2009040)

The reporting follows the recommendations of:

> Charniga K, et al. *Best practices for estimating and reporting
> epidemiological delay distributions of infectious diseases.* 2024.
> [arXiv:2405.08841](https://arxiv.org/abs/2405.08841)

## Authors

Sebastian Funk, Sam Abbott (London School of Hygiene & Tropical Medicine).

## License

MIT (see [LICENSE](LICENSE)).
