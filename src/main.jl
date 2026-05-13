## Package entry points.
## - `fit_joint(d, edges; ...)` runs NUTS on `joint_model`, returns the chain.
##   Shared by `analyse()` and any external code that wants a fit without
##   the full IO pipeline.
## - `analyse(; ...)` loads data, fits, summarises, and writes outputs.
## - `(@main)(args)` is the CLI entry point invoked by `julia -m Hantavirus`.

"""
    fit_joint(d, edges;
              samples = 1000,
              chains  = max(2, min(Threads.nthreads(), 4)),
              seed    = nothing,
              progress = true,
              foffspring_alg = _F_OFFSPRING_ALG)

Run NUTS on `joint_model(d, edges, foffspring_alg)` with Mooncake
reverse-mode AD and return the `MCMCChains.Chains` object. The sampler
configuration (AD backend, NUTS target, init scheme, threading) is
fixed here so `analyse()` and the documentation walk-through share one
source of truth.
"""
function fit_joint(d, edges;
                   samples::Integer = 1000,
                   chains::Integer  = max(2, min(Threads.nthreads(), 4)),
                   seed::Union{Nothing,Integer} = nothing,
                   progress::Bool = true,
                   foffspring_alg = _F_OFFSPRING_ALG)
    seed === nothing || Random.seed!(seed)
    # Mooncake reverse-mode AD: consumes Integrals.jl's ChainRules rrule
    # for `__solvebp`, which lets F_offspring sit on the gradient path
    # without an EnzymeRule.
    adtype = AutoMooncake(; config = Mooncake.Config())
    return sample(
        joint_model(d, edges, foffspring_alg),
        NUTS(0.95; adtype), MCMCThreads(), samples, chains;
        initial_params = fill(DynamicPPL.InitFromPrior(), chains),
        progress = progress,
    )
end

function analyse(;
    data     = LINELIST_PATH,
    obs_time::Union{Nothing,Date} = nothing,
    t0::Union{Nothing,Date} = nothing,
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

    chn = fit_joint(d, edges; samples, chains, seed, progress, foffspring_alg)

    post = summarise(chn)
    save_posterior(post, joinpath(output, "posterior.csv"))
    plot_rt(post, joinpath(figures, "Rt.png"))
    plot_delta_sense_check(chn, d, joinpath(figures, "delta_sense_check.png"))
    plot_pairplot(post, joinpath(figures, "pairplot.png"))
    plot_prior_predictives(joinpath(figures, "prior_predictives.png"))
    plot_posterior_predictions(chn, d, joinpath(figures, "posterior_predictions.png"))
    return chn, post
end

function (@main)(args)
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
    return 0  # exit code expected by `julia -m`
end
