# Getting started

This page walks through a minimal end-to-end example: loading the line list,
fitting the joint model, and inspecting the R(t) estimate.

## Prerequisites

Julia 1.11 or later is required.
Clone the repository and start Julia in the project root:

```
git clone https://github.com/sbfnk/hantavirus.git
cd hantavirus
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Running from the REPL

Load the package and call `analyse()` with reduced settings for a quick run:

```julia
using Hantavirus

# Quick run: 2 chains, 200 samples each (~1–2 minutes)
analyse(chains = 2, samples = 200, output = "output/", figures = "figures/")
```

The full run (4 chains × 1 000 samples) takes a few minutes on a modern
laptop.

## Running from the command line

```
julia --project=. -t auto -m Hantavirus -- \
    --chains 2 --samples 200 \
    --output output/ --figures figures/
```

## Output

After the run finishes you will find:

- `output/posterior.csv` — thinned posterior draws for all parameters.
- `figures/Rt.png` — time-varying reproduction number R(t) with 95% credible
  interval.
- `figures/pairplot.png` — posterior pairplot of population parameters.
- `figures/gi_si_comparison.png` — empirical GI and SI versus analytical
  predictions.
- `figures/incubation_check.png` — empirical incubation period versus fitted
  LogNormal.

## Inspecting results

The posterior CSV can be read back with any Julia CSV library:

```julia
using CSV, DataFrames

post = CSV.read("output/posterior.csv", DataFrame)

# Median incubation mean and 95% CrI
using Statistics
μ = median(post.mu_inc)
lo, hi = quantile(post.mu_inc, [0.025, 0.975])
println("Incubation mean: $μ d ($lo – $hi)")
```
