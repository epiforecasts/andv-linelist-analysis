# Real-time vs retrospective monitoring

The base `joint_model` fits a closed-out outbreak with complete observation.
In real time, three biases need correcting:

1. Long-incubation cases infected before the cut-off may not yet have
   developed symptoms — observed incubation periods are enriched for
   short delays.
2. Late transmissions from any source may not yet have happened or have
   not yet been linked — observed transmission timings (δ) are enriched
   for early / pre-symptomatic events.
3. Recent source cases have not had time to seed all their offspring —
   the observed offspring count is a downward-biased estimate of R(t)
   near the cut-off.

This page walks through fitting the same outbreak retrospectively and
again at a mid-outbreak cut-off and compares the resulting R(t)
trajectories.

## Setting up the two analyses

The whole real-time machinery is gated on `obs_time` in the data object.
Pass nothing and the model is identical to the retrospective form; pass a
date and the per-case truncation terms and the cluster-completeness
adjustment switch on.

`filter_realtime` simulates the analyst's view at a cut-off: cases with
onset after the cut-off are dropped and the offspring count column is
re-derived so the line list is internally consistent.

```@example realtime
using Hantavirus
using Dates: Date
using Random: Random

obs_date = Date("2018-12-31")

ll      = load_linelist()
ll_rt   = filter_realtime(ll, obs_date)

(retrospective = nrow(ll), realtime = nrow(ll_rt))
```

## The cluster-completeness integral

`F_cluster(Δ, μ_inc, σ_inc, μ_δ, σ_δ)` is the probability that the full
`Inc(src) + δ + Inc(sec)` chain fits in `Δ`. It thins the NB offspring
mean for each source case so that the observed count is calibrated
against the fraction of its true offspring already detectable at the
cut-off.

```@example realtime
using Plots: plot, plot!, hline!, savefig

Δs = 0.0:1.0:120.0
F  = F_cluster.(Δs, 3.0, 0.5, 0.0, 1.0)
plt = plot(Δs, F;
           xlabel = "Δ = obs_time − T_inf(src)  (days)",
           ylabel = "F_cluster",
           legend = false, linewidth = 2,
           title  = "Probability the full chain has completed by Δ")
plt
```

`F_cluster` is a 40-by-40 tensor Gauss-Hermite quadrature in standardised
normal coordinates. It is smooth and statically loopy, which lets Enzyme
differentiate it cleanly through `DifferentiationInterface.jl` — the
gradient evaluation NUTS needs on every step. An `Integrals.jl` /
`HCubatureJL` adaptive reference (`F_cluster_quadrature`) is also
available for unit-testing.

## Fitting both views

The two fits use identical priors; the only difference is whether
`obs_time` is supplied. Small sample budgets are used here so the page
builds quickly; for a publication-quality run, scale up `samples` and
`chains`.

```@example realtime
tmp_out = mktempdir()
tmp_fig = mktempdir()

chn_retro, post_retro = analyse(;
    data    = ll,
    samples = 300, chains = 2, seed = 20260512,
    output  = tmp_out, figures = tmp_fig,
)
nothing # hide
```

```@example realtime
chn_rt, post_rt = analyse(;
    data     = ll,
    obs_time = obs_date,
    samples  = 300, chains = 2, seed = 20260512,
    output   = tmp_out, figures = tmp_fig,
)
nothing # hide
```

## R(t) comparison

The two posteriors should agree on R(t) early in the outbreak (before
the cut-off) but the retrospective fit will recover R(t) for the full
window, while the real-time fit only resolves bins covered by retained
cases. The cluster-completeness correction should pull the real-time
estimate up in the bins immediately before the cut-off, where
retrospective offspring counts are biased downward by clusters whose
chains had not yet completed.

```@example realtime
using Statistics: median

rt_q = post -> (lo  = [quantile(exp.(post.log_R_chain[b]), 0.10) for b in eachindex(post.log_R_chain)],
                med = [quantile(exp.(post.log_R_chain[b]), 0.50) for b in eachindex(post.log_R_chain)],
                hi  = [quantile(exp.(post.log_R_chain[b]), 0.90) for b in eachindex(post.log_R_chain)])

retro = rt_q(post_retro)
rt    = rt_q(post_rt)
bins  = 1:length(retro.med)

plt = plot(; xlabel = "Bin index", ylabel = "R(t) (80% CrI)",
             title  = "Retrospective vs real-time R(t)",
             legend = :topright, ylims = (0.0, 4.0))
plot!(plt, bins, retro.med; ribbon = (retro.med .- retro.lo, retro.hi .- retro.med),
      label = "retrospective", linewidth = 2, marker = :circle, color = :steelblue)
plot!(plt, bins, rt.med;    ribbon = (rt.med    .- rt.lo,    rt.hi    .- rt.med),
      label = "real-time ($(obs_date))", linewidth = 2, marker = :diamond, color = :darkorange)
hline!(plt, [1.0]; linestyle = :dash, color = :grey, label = "")
plt
```

## Caveats

The real-time corrections handle three specific biases — long-incubation
cases, late transmissions, and incomplete clusters. The following are
explicitly **not** corrected:

- geographic or severity surveillance biases,
- the onset-to-report delay (only chain completion is modelled),
- general under-ascertainment,
- incomplete source attribution among the cases that are observed.

In settings where reporting delays dominate, the cluster-completeness
adjustment will under-correct: the model assumes that any onset by the
cut-off is reported by the cut-off.
