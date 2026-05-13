## CLI entry point — called by `julia -m TransmissionLinelist`.
## For interactive use, call analyse() directly with keyword arguments.

"""
    sample_fit(model; samples=1000, chains=4, target_accept=0.95,
               seed=20260508, progress=false)

Run NUTS on `model` using the package's default Enzyme AD backend and
`InitFromPrior()` chain initialisation. Returns the FlexiChain.
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

function analyse(;
    data     = LINELIST_PATH,
    output   = OUTPUT_DIR,
    figures  = FIGURES_DIR,
    samples  = 1000,
    chains   = 4,
    seed     = 20260508,
    progress = true,
)
    ll    = load_linelist(data)
    d     = build_data(ll)
    edges = bin_edges_day(d.t0)
    @info "Loaded line list" n_cases=d.N n_sources=sum(>(0), d.source_idx)

    chn = sample_fit(joint_model(d, edges);
        samples  = samples,
        chains   = chains,
        seed     = seed,
        progress = progress,
    )

    post = summarise(chn)
    save_posterior(post, joinpath(output, "posterior.csv"))

    mkpath(figures)
    _save_makie_figure(plot_rt(chn),                       joinpath(figures, "Rt.png"))
    _save_makie_figure(plot_delta_sense_check(chn, d),     joinpath(figures, "delta_sense_check.png"))
    _save_makie_figure(plot_inc_sense_check(chn, d),       joinpath(figures, "inc_sense_check.png"))
    _save_makie_figure(plot_z_ppc(chn, d),                 joinpath(figures, "z_ppc.png"))
    _save_makie_figure(plot_prior_predictives(),           joinpath(figures, "prior_predictives.png"))
    _save_makie_figure(plot_predictive_distributions(chn), joinpath(figures, "predictive_distributions.png"))
    _save_makie_figure(plot_pair(chn),                     joinpath(figures, "pairplot.png"))

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
