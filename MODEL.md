# Model

The model estimates four quantities jointly from the Epuyén line list: the incubation period of Andes hantavirus, the timing of each documented onward transmission relative to its source case's symptom onset, a weekly time-varying reproduction number R(t), and the over-dispersion of offspring counts per case.

Because exposure and symptom-onset dates in the line list are recorded as windows rather than exact dates, each case is given two latent variables: a true infection time `T_inf` and a true symptom-onset time `T_onset`.
The model jointly explores these latent times together with the parameters governing the four quantities above.

Generation and serial intervals are not given separate priors.
They are computed from the fitted incubation period and transmission timing distributions in post-processing: the generation interval is the transmission timing plus the source's incubation period, and the serial interval is the transmission timing plus the secondary's incubation period.

See the [README](https://github.com/sbfnk/hantavirus/blob/main/README.md) for headline results and [LIMITATIONS.md](https://github.com/sbfnk/hantavirus/blob/main/LIMITATIONS.md) for known caveats.

## Reproduction number

The `R(t)` reported here is the case reproduction number indexed by source symptom onset: the expected number of secondary infections produced by a case whose onset is at time `t`.
Onset is the natural choice for these data: fitted transmission timing is tightly clustered around source onset (`μ_δ ≈ 0`, `σ_δ ≈ 0.6 d`), so a case's offspring are infected within roughly a day of the case becoming symptomatic.
Because transmission is so concentrated at onset, this mostly coincides with the instantaneous reproduction number indexed by infection date.

## Submodels

The joint model decomposes into independent components.
Each subsection below states the math first and names the code entry point at the end.
Notation: ``\mathrm{Normal}_+(\mu, \sigma)`` denotes a Normal truncated to non-negative values.

### Incubation period

```math
\begin{aligned}
\mu_{\mathrm{inc}} &\sim \mathrm{Normal}(3.0,\ 0.5) \\
\sigma_{\mathrm{inc}} &\sim \mathrm{Normal}_+(0,\ 0.5) \\
\mathrm{Inc} &\sim \mathrm{LogNormal}(\mu_{\mathrm{inc}},\ \sigma_{\mathrm{inc}}).
\end{aligned}
```

Implemented in `incubation_model`.

### Transmission timing

The per-pair gap ``\delta = T_{\mathrm{inf}}(\mathrm{sec}) - T_{\mathrm{onset}}(\mathrm{src})`` can be negative (pre-symptomatic transmission).

```math
\begin{aligned}
\mu_\delta &\sim \mathrm{Normal}(0,\ 5) \\
\sigma_\delta &\sim \mathrm{Normal}_+(0,\ 1) \\
\delta &\sim \mathrm{Normal}(\mu_\delta,\ \sigma_\delta).
\end{aligned}
```

Implemented in `transmission_delta_model`.

### Latent infection and onset times

Each case ``i`` has two continuous latents.
``T_{\mathrm{onset}}[i]`` is uniform over the recorded onset window ``[L_i, U_i]``.
``T_{\mathrm{inf}}[i]`` is uniform over the exposure window for sourced cases (or an 80-day pre-onset window for the index), capped above by ``T_{\mathrm{onset}}[i]``:

```math
\begin{aligned}
T_{\mathrm{onset}}[i] &\sim \mathrm{Uniform}(L_i,\ U_i) \\
T_{\mathrm{inf}}[i] &\sim
\begin{cases}
\mathrm{Uniform}\bigl(L_i - 80,\ T_{\mathrm{onset}}[i]\bigr) & \mathrm{src}(i) = 0 \\
\mathrm{Uniform}\bigl(\ell_i,\ \min(u_i,\ T_{\mathrm{onset}}[i])\bigr) & \mathrm{src}(i) \ne 0
\end{cases}
\end{aligned}
```

where ``[\ell_i, u_i]`` is the recorded exposure window.
A per-pair constraint ``T_{\mathrm{inf}}[i] > T_{\mathrm{inf}}(\mathrm{src}(i))`` enforces a positive generation interval.
For each case ``i`` the submodel adds

```math
\log \mathrm{pdf}\bigl(\mathrm{Inc},\ T_{\mathrm{onset}}[i] - T_{\mathrm{inf}}[i]\bigr)
```

and, for sourced cases, also

```math
\log \mathrm{pdf}\bigl(\delta,\ T_{\mathrm{inf}}[i] - T_{\mathrm{onset}}(\mathrm{src}(i))\bigr).
```

Implemented in `latent_times_model`.

### Real-time truncation

Active only when an `obs_time` cut-off is set.
Two contributions.

First, the right-truncation on the index Inc, on the observation event ``T_{\mathrm{inf}}[i] + \mathrm{Inc} \le \mathrm{obs\_time}``:

```math
-\sum_{i: \mathrm{src}(i) = 0} \log \mathrm{cdf}\bigl(\mathrm{Inc},\ \mathrm{obs\_time} - T_{\mathrm{inf}}[i]\bigr).
```

Sourced cases need no such factor; their exposure window already bounds ``T_{\mathrm{inf}}[i]`` below ``T_{\mathrm{onset}}[i] \le \mathrm{obs\_time}``.

Second, the per-pair offspring-completeness denominator.
With ``p_j = \mathrm{cdf}\bigl(\delta + \mathrm{Inc}(\mathrm{sec}),\ \mathrm{obs\_time} - T_{\mathrm{onset}}[j]\bigr)``, each sourced case ``i`` contributes

```math
-\log p_{\mathrm{src}(i)}.
```

Implemented in `truncation_model`; ``p`` uses the `ConvolvedDelays` distribution.

### Reproduction number

``\log R`` evolves as a random walk on weekly knots ``b = 1, \ldots, B``:

```math
\begin{aligned}
\log R_1 &\sim \mathrm{Normal}(\log 1.5,\ 1) \\
\sigma_{\mathrm{rw}} &\sim \mathrm{Normal}_+(0,\ 0.5) \\
\varepsilon_b &\sim \mathrm{Normal}(0,\ 1), \quad b = 1, \ldots, B - 1 \\
\log R_b &= \log R_1 + \sigma_{\mathrm{rw}} \sum_{j=1}^{b-1} \varepsilon_j.
\end{aligned}
```

Between knots, ``\log R(t)`` is linearly interpolated.
The reported ``R(t)`` is the *case reproduction number* indexed by source symptom onset.
Implemented in `random_walk_rt_model`.

### Negative-Binomial dispersion

The standard ``1/\sqrt k`` reparameterisation:

```math
\begin{aligned}
\phi^{-1/2} &\sim \mathrm{Normal}_+(0,\ 1) \\
k &= 1 / (\phi^{-1/2})^2.
\end{aligned}
```

Implemented in `nb_dispersion_model`.

### Case likelihood

For each case ``i`` with onset ``T_{\mathrm{onset}}[i]`` falling in bin ``b(i)``:

```math
Z_{\mathrm{obs}}[i] \sim \mathrm{NegBin}\!\bigl(k,\ \exp(\log R_{b(i)})\cdot p_i\bigr)
```

with the offspring-completeness ``p_i = 1`` in retrospective mode and ``p_i = \mathrm{cdf}\bigl(\delta + \mathrm{Inc}(\mathrm{sec}),\ \mathrm{obs\_time} - T_{\mathrm{onset}}[i]\bigr)`` in real-time mode.
The Negative-Binomial is mean–dispersion parameterised so ``\mathbb{E}[Z_{\mathrm{obs}}[i]] = \exp(\log R_{b(i)})\cdot p_i`` and ``\mathrm{Var}[Z_{\mathrm{obs}}[i]] = \exp(\log R_{b(i)})\cdot p_i \cdot (1 + \exp(\log R_{b(i)})\cdot p_i / k)``.
Implemented in `case_model`.

## Inference

NUTS, 4 chains, 1000 post-warmup samples each, `target_accept = 0.95`.
Default seed: 20260508.
Reverse-mode AD via Mooncake.

## Real-time predictions

Two counterfactual predictors give posterior-predictive future onset counts conditional on a real-time fit at cut-off `obs_time`.
Both reuse the same Gamma–Poisson conjugate update on each source's true offspring rate.

For source `i` with onset at `T_onset[i]` let

```math
\Delta_i = \mathrm{obs\_time} - T_{\mathrm{onset}}[i],\quad
p_i = \Pr(\delta + \mathrm{Inc}(\mathrm{sec}) \le \Delta_i),\quad
q_i = \Pr(\delta \le \Delta_i),
```

where ``p_i`` is the offspring-completeness probability (chain finished by `obs_time`, the cdf of `ConvolvedDelays(inc, δ)`) and ``q_i`` is the probability that transmission has happened by `obs_time` (the cdf of the transmission timing distribution).

The latent rate of source ``i``'s offspring follows the conjugate posterior

```math
\lambda_i \mid Z_{\mathrm{obs}}[i], k, R_i, p_i \;\sim\; \mathrm{Gamma}\!\left(k + Z_{\mathrm{obs}}[i],\ \frac{R_i}{k + R_i\,p_i}\right) \quad (\text{scale form}).
```

The two predictors differ only in the thinning probability applied to ``\lambda_i`` when drawing future onsets:

- **Controlled** (`predict_controlled_outbreak`): transmission stops at `obs_time`.
  Only people already infected by then (`δ ≤ Δ_i`, chain not yet symptomatic) contribute:

```math
Z_{\mathrm{future}}[i] \;\sim\; \mathrm{Poisson}\!\bigl(\lambda_i\,(q_i - p_i)\bigr).
```

- **Natural chain** (`predict_natural_chain_outbreak`): current sources keep transmitting at their existing rate but no second-generation chains form from those new offspring:

```math
Z_{\mathrm{future}}[i] \;\sim\; \mathrm{Poisson}\!\bigl(\lambda_i\,(1 - p_i)\bigr).
```

Both sum the contributions across observed sources for each posterior draw.
The strict counterfactual is a subset of the natural-chain one: every offspring in the controlled prediction is also in the natural-chain prediction, plus the offspring whose transmission hasn't happened yet.
