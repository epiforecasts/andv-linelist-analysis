module Hantavirus

using ArgParse: ArgParseSettings, @add_arg_table!, parse_args
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
include("main.jl")

export load_linelist, build_data, bin_edges_day, which_bin, bin_labels
export LINELIST_PATH, OUTPUT_DIR, BIN_EDGES
export joint_model
export diagnostics, vector_chain, summarise, save_posterior
export main

end
