module Hantavirus

using CSV: CSV
using DataFrames: DataFrame, nrow, eachrow, passmissing
using Dates: Dates, Date, Day
using Distributions: Normal, LogNormal, truncated, NegativeBinomial,
                     Uniform, logpdf, cdf
using MCMCChains: MCMCChains
using Printf: @printf, @sprintf
using Random: Random
using Statistics: quantile
using Turing: Turing, @model, NUTS, MCMCThreads, sample
import FlexiChains

include("data.jl")
include("model.jl")
include("postprocess.jl")

function main(args = ARGS; n_chains = 4, n_samples = 1000, seed = 20260508)
    Random.seed!(seed)

    ll    = load_linelist()
    data  = build_data(ll)
    edges = bin_edges_day(data.t0)
    @info "Loaded line list" n_cases=data.N n_sources=sum(>(0), data.source_idx)

    chn = sample(
        joint_model(data, edges),
        NUTS(0.95), MCMCThreads(), n_samples, n_chains;
        progress = false,
    )

    post = summarise(chn)
    save_posterior(post, joinpath(OUTPUT_DIR, "posterior.csv"))
    return chn, post
end

export load_linelist, build_data, bin_edges_day, which_bin, bin_labels
export LINELIST_PATH, OUTPUT_DIR, BIN_EDGES
export joint_model
export diagnostics, vector_chain, summarise, save_posterior
export main

end
