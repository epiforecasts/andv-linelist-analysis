## CLI entry point — called by `julia -m Hantavirus`.

using ArgParse: ArgParseSettings, @add_arg_table!, parse_args

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

    Random.seed!(p["seed"])

    ll    = load_linelist(p["data"])
    data  = build_data(ll)
    edges = bin_edges_day(data.t0)
    @info "Loaded line list" n_cases=data.N n_sources=sum(>(0), data.source_idx)

    chn = sample(
        joint_model(data, edges),
        NUTS(0.95), MCMCThreads(), p["samples"], p["chains"];
        progress = false,
    )

    post = summarise(chn)
    save_posterior(post, joinpath(p["output"], "posterior.csv"))
    return chn, post
end
