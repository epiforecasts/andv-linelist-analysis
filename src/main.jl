## Package entry points.
## - `sample_fit(model; ...)` runs NUTS on any Turing model.
## - `analyse(; ...)` loads data, fits, summarises, and writes outputs.
## - `main(args)` is the CLI entry point invoked by `julia -m Hantavirus`.

"""
    sample_fit(model; samples = 1000, chains = ..., target_accept = 0.95,
               seed = nothing, progress = false,
               adtype = AutoMooncake(; config = Mooncake.Config()))

Run NUTS on `model` and return the `MCMCChains.Chains` object.

Mooncake reverse-mode AD is the default because it is the only backend
that survives `Integrals.jl`'s `solve` through the real-time
`F_offspring` cubature; Enzyme works for the retrospective form but
trips inside HCubature for the real-time form.
"""
function sample_fit(model;
    samples::Integer       = 1000,
    chains::Integer        = max(2, min(Threads.nthreads(), 4)),
    target_accept::Real    = 0.95,
    seed::Union{Nothing,Integer} = nothing,
    progress::Bool         = false,
    adtype = AutoMooncake(; config = Mooncake.Config()),
)
    seed === nothing || Random.seed!(seed)
    return sample(
        model,
        NUTS(target_accept; adtype), MCMCThreads(), samples, chains;
        initial_params = fill(DynamicPPL.InitFromPrior(), chains),
        progress = progress,
    )
end

function analyse(;
    data     = LINELIST_PATH,
    obs_time::Union{Nothing,Date} = nothing,
    t0::Union{Nothing,Date}       = nothing,
    output   = OUTPUT_DIR,
    figures  = FIGURES_DIR,
    samples  = 1000,
    chains   = max(2, min(Threads.nthreads(), 4)),
    seed     = 20260508,
    progress = true,
    foffspring_alg = _F_OFFSPRING_ALG,
)
    ll = data isa DataFrame ? data : load_linelist(data)
    if obs_time !== nothing
        ll = filter_realtime(ll, obs_time)
    end
    d     = build_data(ll; obs_time = obs_time, t0 = t0)
    edges = bin_edges_day(d.t0)
    @info "Loaded line list" n_cases=d.N n_sources=sum(>(0), d.source_idx) obs_time=obs_time

    chn = sample_fit(joint_model(d, edges, foffspring_alg);
                     samples, chains, seed, progress)

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

function main(args)
    s = ArgParseSettings(; description = "Fit joint ANDV incubation/R(t) model")
    @add_arg_table! s begin
        "--data", "-d"
            help    = "path to linelist CSV"
            default = LINELIST_PATH
        "--output", "-o"
            help    = "output directory for posterior.csv"
            default = OUTPUT_DIR
        "--figures", "-f"
            help    = "directory for figures"
            default = FIGURES_DIR
        "--samples", "-n"
            help     = "NUTS samples per chain"
            arg_type = Int
            default  = 1000
        "--chains", "-c"
            help     = "number of parallel chains (default: clamp(Threads.nthreads(), 2, 4))"
            arg_type = Int
            default  = max(2, min(Threads.nthreads(), 4))
        "--seed", "-s"
            help     = "random seed"
            arg_type = Int
            default  = 20260508
        "--obs-time"
            help    = "real-time cut-off date (ISO format, e.g. 2018-12-31); " *
                      "omit for a retrospective fit"
            default = nothing
    end
    p = parse_args(args, s)
    obs_time = p["obs-time"] === nothing ? nothing : Date(p["obs-time"])
    analyse(;
        data     = p["data"],
        obs_time = obs_time,
        output   = p["output"],
        figures  = p["figures"],
        samples  = p["samples"],
        chains   = p["chains"],
        seed     = p["seed"],
        progress = false,
    )
    return 0
end
