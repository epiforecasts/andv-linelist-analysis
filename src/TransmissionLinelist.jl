module TransmissionLinelist

using ArgParse: ArgParseSettings, @add_arg_table!, parse_args
using Chain: @chain
using CSV: CSV
using DataFrames: DataFrame, nrow, eachrow, passmissing, rename!, sort!
using DataFramesMeta: @select, @transform, @subset, @combine, @by, @rtransform,
                     @rsubset, @orderby, @rename
using Dates: Dates, Date, Day, dayofweek
using Distributions: Normal, LogNormal, truncated, NegativeBinomial,
                     Uniform, logpdf, cdf, pdf
using DocStringExtensions: TYPEDSIGNATURES
using PairPlots: PairPlots, pairplot
using Makie: Makie, Figure, Axis, Legend, Theme, Auto, BarPlot,
             theme_latexfonts, with_theme,
             lines!, scatter!, hist!, hlines!, vlines!, band!,
             axislegend, rowsize!
using AlgebraOfGraphics: AlgebraOfGraphics, data, mapping, visual, draw!
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
export prepare_model
export joint_model
export diagnostics, diagnostics_table, summary_table
export vector_chain, summarise, save_posterior
export plot_data, plot_rt, plot_pair, plot_posterior_predictive
export plot_delta_sense_check, plot_prior_predictives
export analyse, sample_fit

end
