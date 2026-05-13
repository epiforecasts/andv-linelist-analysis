## CLI entry point — called by `julia -m TransmissionLinelist`.
## For interactive use, call analyse() directly with keyword arguments.

"""
$(TYPEDSIGNATURES)

Run NUTS on `model` using Enzyme AD and `InitFromPrior()` initialisation.

Wraps `Turing.sample` with the package's defaults: Enzyme reverse-mode
gradients, multi-chain `MCMCThreads()` and prior-based chain
initialisation. Used by [`analyse`](@ref) and the analysis walkthrough.

# Arguments
- `model`: a `Turing` model, e.g. the output of [`joint_model`](@ref).

# Keyword Arguments
- `samples`: number of NUTS samples per chain. Defaults to `1000`.
- `chains`: number of parallel chains. Defaults to `4`.
- `target_accept`: NUTS target acceptance rate. Defaults to `0.95`.
- `seed`: random seed passed to `Random.seed!`. Defaults to `20260508`.
- `progress`: whether to display the sampler's progress bar. Defaults
  to `false`.

# Returns
A FlexiChain (as returned by `Turing.sample`).
"""
function sample_fit(model;
    samples       = 1000,
    chains        = 4,
    target_accept = 0.95,
    seed          = 20260508,
    progress      = false,
)
    Random.seed!(seed)
    adtype = AutoEnzyme(; mode = Enzyme.set_runtime_activity(Enzyme.Reverse))
    return sample(
        model,
        NUTS(target_accept; adtype), MCMCThreads(), samples, chains;
        initial_params = fill(DynamicPPL.InitFromPrior(), chains),
        progress = progress,
    )
end

"""
$(TYPEDSIGNATURES)

Fit the joint ANDV model and write posterior summaries and figures.

Loads the line list, builds the model data, runs NUTS with Enzyme
gradients across multiple chains, prints the [`summarise`](@ref) output,
writes `posterior.csv` to `output`, and renders five PNG figures into
`figures` (`Rt.png`, `delta_sense_check.png`, `pairplot.png`,
`prior_predictives.png`, `posterior_predictions.png`). This is the
single source of truth for the CLI options exposed by [`main`](@ref)
and the Options table in the README.

# Keyword Arguments
- `data`: path to the line-list CSV. Defaults to the bundled
  `data/linelist.csv`.
- `output`: directory in which to write `posterior.csv`. Defaults to
  the package's `output/` directory.
- `figures`: directory in which to write the figures. Defaults to the
  package's `figures/` directory.
- `samples`: number of NUTS samples per chain. Defaults to `1000`.
- `chains`: number of parallel chains. Defaults to `4`.
- `seed`: random seed. Defaults to `20260508`.
- `progress`: whether to display the sampler's progress bar. Defaults
  to `true`.

# Returns
A `(chn, post)` tuple containing the raw `Turing.sample` chain and the
summary tuple returned by [`summarise`](@ref).

# Examples
```julia
using TransmissionLinelist
chn, post = analyse(samples = 500, chains = 2)
```
"""
function analyse(;
    data     = LINELIST_PATH,
    output   = OUTPUT_DIR,
    figures  = FIGURES_DIR,
    samples  = 1000,
    chains   = 4,
    seed     = 20260508,
    progress = true,
)
    ll = load_linelist(data)
    model, d, _ = prepare_model(ll)
    @info "Loaded line list" n_cases=d.N n_sources=sum(>(0), d.source_idx)

    chn = sample_fit(model;
        samples  = samples,
        chains   = chains,
        seed     = seed,
        progress = progress,
    )

    post = summarise(chn)
    save_posterior(post, joinpath(output, "posterior.csv"))

    mkpath(figures)
    _save_makie_figure(plot_rt(chn),                   joinpath(figures, "Rt.png"))
    _save_makie_figure(plot_delta_sense_check(chn, d), joinpath(figures, "delta_sense_check.png"))
    _save_makie_figure(plot_prior_predictives(),       joinpath(figures, "prior_predictives.png"))
    _save_makie_figure(plot_posterior_predictive(chn), joinpath(figures, "posterior_predictions.png"))
    _save_makie_figure(plot_pair(chn),                 joinpath(figures, "pairplot.png"))

    return chn, post
end

# Saving a Makie figure needs a Makie backend (e.g. CairoMakie) loaded at
# the call site. Look the `save` method up dynamically so the package
# itself doesn't depend on a particular backend.
function _save_makie_figure(fig, path)
    backend = nothing
    for name in (:CairoMakie, :GLMakie, :WGLMakie)
        if isdefined(Main, name)
            backend = getfield(Main, name)
            break
        end
    end
    if backend === nothing
        @warn "No Makie backend loaded in Main; skipping figure save" path
        return path
    end
    Base.invokelatest(backend.save, path, fig; px_per_unit = 2.0)
    return path
end

"""
$(TYPEDSIGNATURES)

Command-line entry point invoked by `julia -m TransmissionLinelist`.

Parses `args` (typically `ARGS`) with `ArgParse` and forwards the parsed
options to [`analyse`](@ref). The `help=` strings here mirror the
`# Keyword Arguments` section of [`analyse`](@ref), which is the single
source of truth.

# Arguments
- `args`: a vector of CLI argument strings (typically `ARGS`).

# Returns
`0` as the process exit code, as expected by Julia's `julia -m` entry
point convention.
"""
function (@main)(args)
    s = ArgParseSettings(; description = "Fit joint ANDV incubation/R(t) model")
    @add_arg_table! s begin
        "--data", "-d"
            help    = "path to the line-list CSV"
            default = LINELIST_PATH
        "--output", "-o"
            help    = "directory in which to write posterior.csv"
            default = OUTPUT_DIR
        "--figures", "-f"
            help    = "directory in which to write the figures"
            default = FIGURES_DIR
        "--samples", "-n"
            help     = "number of NUTS samples per chain"
            arg_type = Int
            default  = 1000
        "--chains", "-c"
            help     = "number of parallel chains"
            arg_type = Int
            default  = 4
        "--seed", "-s"
            help     = "random seed"
            arg_type = Int
            default  = 20260508
    end
    p = parse_args(args, s)
    analyse(;
        data     = p["data"],
        output   = p["output"],
        figures  = p["figures"],
        samples  = p["samples"],
        chains   = p["chains"],
        seed     = p["seed"],
        progress = false,
    )
    return 0  # exit code expected by `julia -m`
end
