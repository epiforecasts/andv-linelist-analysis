## CLI entry point — called by `julia -m Hantavirus`.
## For interactive use, call analyse() directly with keyword arguments.

function analyse(;
    data     = LINELIST_PATH,
    output   = OUTPUT_DIR,
    figures  = FIGURES_DIR,
    samples  = 1000,
    chains   = 4,
    seed     = 20260508,
    progress = true,
)
    Random.seed!(seed)

    ll    = load_linelist(data)
    d     = build_data(ll)
    edges = bin_edges_day(d.t0)
    @info "Loaded line list" n_cases=d.N n_sources=sum(>(0), d.source_idx)

    adtype = AutoEnzyme(; mode = Enzyme.set_runtime_activity(Enzyme.Reverse))
    chn = sample(
        joint_model(d, edges),
        NUTS(0.95; adtype), MCMCThreads(), samples, chains;
        initial_params = fill(DynamicPPL.InitFromPrior(), chains),
        progress = progress,
    )

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
