## CLI entry point — called by `julia -m Hantavirus`.
## For interactive use, call analyse() directly with keyword arguments.

function analyse(;
    data     = LINELIST_PATH,
    output   = OUTPUT_DIR,
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

    chn = sample(
        joint_model(d, edges),
        NUTS(0.95), MCMCThreads(), samples, chains;
        progress = progress,
    )

    post = summarise(chn)
    save_posterior(post, joinpath(output, "posterior.csv"))
    figures = joinpath(output, "figures")
    plot_rt(post, joinpath(figures, "Rt.png"))
    plot_delta_sense_check(chn, d, joinpath(figures, "delta_sense_check.png"))
    return chn, post
end

function main(args = ARGS)
    s = ArgParseSettings(; description = "Fit joint ANDV incubation/R(t) model")
    @add_arg_table! s begin
        "--data", "-d"
            help    = "path to linelist CSV"
            default = LINELIST_PATH
        "--output", "-o"
            help    = "output directory for posterior.csv"
            default = OUTPUT_DIR
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
    return analyse(;
        data     = p["data"],
        output   = p["output"],
        samples  = p["samples"],
        chains   = p["chains"],
        seed     = p["seed"],
        progress = false,
    )
end
