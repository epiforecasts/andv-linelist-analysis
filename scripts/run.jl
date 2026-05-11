## Entry point: load the line list, fit the joint Inc + δ + R(t) model,
## summarise, and save the posterior to output/posterior.csv.
##
## Run with:  julia --project=. -t auto scripts/run.jl

using Random, Turing

const SRC_DIR = normpath(joinpath(@__DIR__, "..", "src"))
include(joinpath(SRC_DIR, "data.jl"))
include(joinpath(SRC_DIR, "model.jl"))
include(joinpath(SRC_DIR, "postprocess.jl"))

function main(; n_chains = 4, n_samples = 1000, seed = 20260508)
    Random.seed!(seed)

    ll    = load_linelist()
    data  = build_data(ll)
    edges = bin_edges_day(data.t0)
    @info "Loaded line list" n_cases = data.N n_sources = sum(>(0), data.source_idx)

    chn = sample(joint_model(data, edges),
                 NUTS(0.95), MCMCThreads(), n_samples, n_chains; progress = false)

    post = summarise(chn)
    save_posterior(post, joinpath(OUTPUT_DIR, "posterior.csv"))
    plot_rt(post, joinpath(PROJECT_ROOT, "figures", "Rt.png"))
    return chn, post
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
