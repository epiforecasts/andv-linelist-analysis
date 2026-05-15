## Package entry points.
## - `sample_fit(model; ...)` runs NUTS on any Turing model.
## - `analyse(; ...)` loads data, fits, summarises, and writes outputs.
## - `main(args)` is the CLI entry point invoked by `julia -m TransmissionLinelist`.

"""
$(TYPEDSIGNATURES)

Run NUTS on `model` and return the sampled chain.

Uses Mooncake as the default AD backend and `InitFromPrior()` chain
initialisation. Override `adtype` to swap backends (e.g. Enzyme).

# Arguments
- `model`: a Turing model, e.g. from [`joint_model`](@ref).

# Keyword Arguments
- `samples`: NUTS samples per chain.
- `chains`: number of parallel chains.
- `target_accept`: NUTS acceptance target.
- `seed`: optional random seed.
- `progress`: show a NUTS progress bar.
- `adtype`: AD backend for NUTS. Defaults to
  `AutoMooncake(; config = Mooncake.Config())`.
"""
function sample_fit(model;
        samples::Integer = 1000,
        chains::Integer = 4,
        target_accept::Real = 0.95,
        seed::Integer = 20260508,
        progress::Bool = false,
        adtype = AutoMooncake(; config = Mooncake.Config())
)
    Random.seed!(seed)
    return sample(
        model,
        NUTS(target_accept; adtype), MCMCThreads(), samples, chains;
        initial_params = fill(DynamicPPL.InitFromPrior(), chains),
        progress = progress
    )
end

"""
$(TYPEDSIGNATURES)

Load the line list, fit the joint model, save the posterior summary to
`output/posterior.csv`, and (unless `plots = false`) write all figures
into `figures/`. Returns `(chain, post)`.

# Keyword Arguments
- `data`: path to the line-list CSV, or a pre-loaded `DataFrame`.
- `obs_time`: optional real-time cut-off `Date`; omit for a retrospective fit.
- `t0`: optional explicit time origin (`Date`); defaults to
  `minimum(onset_date) - 60 d`.
- `output`: directory for `posterior.csv`.
- `figures`: directory for the figure PNGs.
- `samples`: NUTS samples per chain.
- `chains`: number of parallel chains.
- `seed`: random seed.
- `progress`: show a NUTS progress bar.
- `plots`: skip all figure generation when `false`.
- `foffspring_alg`: integration algorithm used for `F_offspring` in
  real-time mode.
"""
function analyse(;
        data = LINELIST_PATH,
        obs_time::Union{Nothing, Date} = nothing,
        t0::Union{Nothing, Date} = nothing,
        output = OUTPUT_DIR,
        figures = FIGURES_DIR,
        samples = 1000,
        chains = 4,
        seed = 20260508,
        progress = true,
        plots = true,
        foffspring_alg = _F_OFFSPRING_ALG
)
    ll = data isa DataFrame ? data : load_linelist(data)
    if obs_time !== nothing
        ll = filter_realtime(ll, obs_time)
    end
    d = build_data(ll; obs_time = obs_time, t0 = t0)
    edges = bin_edges_day(d.t0)
    if obs_time !== nothing
        obs_offset = Float64(Dates.value(obs_time - d.t0))
        edges = edges[edges .<= obs_offset]
        if isempty(edges) || edges[end] < obs_offset
            push!(edges, obs_offset)
        end
    end
    @info "Loaded line list" n_cases=d.N n_sources=sum(>(0), d.source_idx) obs_time=obs_time n_knots=length(edges)

    chn = sample_fit(joint_model(d, edges, foffspring_alg);
        samples, chains, seed, progress)

    post = summarise(chn)
    save_posterior(post, joinpath(output, "posterior.csv"))

    if !plots && figures != FIGURES_DIR
        @warn "plots=false; --figures path ignored" figures
    end
    if plots
        mkpath(figures)
        _save_makie_figure(plot_rt(chn), joinpath(figures, "Rt.png"))
        _save_makie_figure(plot_delta_sense_check(chn, d),
            joinpath(figures, "delta_sense_check.png"))
        _save_makie_figure(plot_inc_sense_check(chn, d),
            joinpath(figures, "inc_sense_check.png"))
        _save_makie_figure(plot_z_ppc(chn, d), joinpath(figures, "z_ppc.png"))
        _save_makie_figure(plot_prior_predictives(), joinpath(figures, "prior_predictives.png"))
        _save_makie_figure(plot_predictive_distributions(chn),
            joinpath(figures, "predictive_distributions.png"))
        _save_makie_figure(plot_pair(chn), joinpath(figures, "pairplot.png"))
    end

    return chn, post
end

# Saving a Makie figure needs a Makie backend (e.g. CairoMakie) loaded
# at the call site. Look up `save` dynamically so the package itself
# doesn't depend on a particular backend.
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

CLI entry point invoked by `julia -m TransmissionLinelist`. Parses
command-line arguments and forwards them to [`analyse`](@ref).

# Arguments
- `args`: vector of command-line argument strings.
"""
function main(args)
    s = ArgParseSettings(; description = "Fit joint ANDV incubation/R(t) model")
    @add_arg_table! s begin
        "--data", "-d"
        help = "path to linelist CSV"
        default = LINELIST_PATH
        "--output", "-o"
        help = "output directory for posterior.csv"
        default = OUTPUT_DIR
        "--figures", "-f"
        help = "directory for figures"
        default = FIGURES_DIR
        "--no-figures"
        help = "skip all figure generation"
        action = :store_true
        "--samples", "-n"
        help = "NUTS samples per chain"
        arg_type = Int
        default = 1000
        "--chains", "-c"
        help = "number of parallel chains"
        arg_type = Int
        default = 4
        "--seed", "-s"
        help = "random seed"
        arg_type = Int
        default = 20260508
        "--obs-time"
        help = "real-time cut-off date (ISO format, e.g. 2018-12-31); " *
               "omit for a retrospective fit"
        default = nothing
    end
    p = parse_args(args, s)
    obs_time = p["obs-time"] === nothing ? nothing : Date(p["obs-time"])
    analyse(;
        data = p["data"],
        obs_time = obs_time,
        output = p["output"],
        figures = p["figures"],
        samples = p["samples"],
        chains = p["chains"],
        seed = p["seed"],
        progress = false,
        plots = !p["no-figures"]
    )
    return 0
end
