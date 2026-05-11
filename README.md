# Andes virus — joint estimation of incubation, transmission timing, and R(t)

A Julia + Turing model for the Epuyén 2018–19 Andes hantavirus outbreak
([Martínez et al. 2020, NEJM](https://doi.org/10.1056/NEJMoa2009040)).
From the line list given in the paper it jointly estimates the incubation
period and the transmission timing of each secondary infection relative to
its source's symptom onset, plus a time-varying reproduction number with
offspring dispersion. Double interval censoring is handled by a continuous
latent infection time for each case.

The generation interval (transmission timing plus the source's incubation
period) and the serial interval (transmission timing plus the secondary's
incubation period) are derived in post-processing from the fitted
distributions. A per-pair constraint that the secondary's infection time is
later than the source's keeps the generation interval positive at the
latent level.

## Headline results (Epuyén line list)

### Incubation period (estimated)

| Quantity | Posterior median (95% CrI) |
|---|---|
| Mean | 22.5 d (20.3 – 25.3) |
| 95th percentile | 36.0 d (31.3 – 43.8) |
| 99th percentile | 44.7 d (37.5 – 57.6) |

### Transmission timing relative to source onset (estimated)

Negative values mean the secondary was infected before the source became symptomatic.

| Quantity | Posterior median (95% CrI) |
|---|---|
| Mean | 0.17 d (−0.17 – 0.49) |
| SD | 0.61 d (0.46 – 0.83) |
| P(transmission is pre-symptomatic) | 0.39 (0.20 – 0.61) |
| P(pre-symptomatic by more than 1 day) | 0.03 (0.00 – 0.12) |
| P(pre-symptomatic by more than 2 days) | 0.00 (0.00 – 0.01) |

### Generation interval / serial interval (derived from incubation and transmission timing)

| Quantity | Posterior median (95% CrI) |
|---|---|
| Mean | 22.7 d (20.4 – 25.5) |
| SD | 7.3 d (5.6 – 10.3) |

### Offspring and time-varying reproduction number

| Quantity | Posterior median (95% CrI) |
|---|---|
| Negative-Binomial offspring dispersion | 0.33 (0.12 – 0.88) |
| R(t) — Nov 2018 | 1.32 (0.60 – 3.13) |
| R(t) — Dec 2018 | 0.51 (0.14 – 1.67) |
| R(t) — Jan 2019 | 0.40 (0.04 – 2.04) |
| R(t) — Feb 2019+ | 0.42 (0.02 – 2.99) |

## Repository layout

```
src/
  Hantavirus.jl    — module entry point and imports
  data.jl          — line list loading and bin definitions
  model.jl         — the joint Turing model (incubation, transmission timing, R(t))
  postprocess.jl   — diagnostics, summaries, CSV output
  main.jl          — CLI entry point (argument parsing)
data/
  linelist.csv     — Epuyén outbreak line list (Martínez Table S2)
Project.toml       — Julia package manifest
Manifest.toml      — locked dependency versions
LICENSE            — MIT
```

## Data

The Epuyén line list (`data/linelist.csv`) is hand-encoded from Table S2 of
the supplementary appendix of Martínez et al. 2020.

## Running

```
julia --project=. -t auto -m Hantavirus
```

NUTS, 4 chains × 1000 samples. Takes a few minutes on a laptop. Posterior
saved to `output/posterior.csv`.

Options:

```
-d, --data      path to linelist CSV   (default: data/linelist.csv)
-o, --output    output directory        (default: output/)
-n, --samples   NUTS samples per chain  (default: 1000)
-c, --chains    number of chains        (default: 4)
-s, --seed      random seed             (default: 20260508)
```

Example:

```
julia --project=. -t auto -m Hantavirus -- -n 500 -c 2 -o results/
```

### From the REPL

```julia
julia> using Hantavirus
julia> main()                                          # all defaults
julia> main(["-n", "500", "-c", "2", "-o", "results/"])  # with options
```

## Citing

If you use this code or the Epuyén line list encoding, please cite:

> Martínez VP, Di Paola N, Alonso DO, et al. *"Super-spreaders" and
> person-to-person transmission of Andes virus in Argentina.* N Engl J Med
> 2020;383:2230–41. [doi:10.1056/NEJMoa2009040](https://doi.org/10.1056/NEJMoa2009040)

## License

MIT (see [LICENSE](LICENSE)).
