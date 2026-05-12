module Hantavirus

using ArgParse: ArgParseSettings, @add_arg_table!, parse_args
using Chain: @chain
using CSV: CSV
using DataFrames: DataFrame, nrow, eachrow, passmissing, rename!, sort!
using DataFramesMeta: @select, @transform, @subset, @combine, @by, @rtransform,
                     @rsubset, @orderby, @rename
using Dates: Dates, Date, Day, dayofweek
using Distributions: Normal, LogNormal, truncated, NegativeBinomial,
                     Uniform, logpdf, cdf, pdf
using PairPlots: PairPlots, pairplot
using Plots: plot, plot!, bar, hline!, histogram, histogram!, vline!,
             scatter, scatter!, savefig
using Printf: @printf, @sprintf
using Random: Random
using Statistics: quantile, mean, std
using Turing: Turing, @model, NUTS, MCMCThreads, sample, DynamicPPL
using ADTypes: AutoEnzyme
using Enzyme: Enzyme
import FlexiChains

include("data.jl")
include("model.jl")
include("postprocess.jl")
include("plots.jl")
include("main.jl")

export load_linelist, build_data, bin_edges_day, log_R_at, bin_labels
export joint_model
export diagnostics, diagnostics_table, summary_table
export vector_chain, summarise, save_posterior
export plot_data, plot_rt, plot_pair, plot_posterior_predictive
export plot_delta_sense_check, plot_prior_predictives
export analyse, main

end
