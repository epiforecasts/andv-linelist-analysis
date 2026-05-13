module Hantavirus

using ADTypes: AutoMooncake, AutoForwardDiff
using AlgebraOfGraphics: AlgebraOfGraphics, data, mapping, visual, draw!
using ArgParse: ArgParseSettings, @add_arg_table!, parse_args
using Chain: @chain
using CSV: CSV
using DataFrames: DataFrame, nrow, eachrow, passmissing, rename!, sort!
using DataFramesMeta: @select, @transform, @subset, @combine, @by, @rtransform,
                     @rsubset, @orderby, @rename
using Dates: Dates, Date, Day, dayofweek
using Distributions: Normal, LogNormal, truncated, NegativeBinomial,
                     Uniform, logpdf, logcdf, cdf, pdf
using Integrals: IntegralProblem, GaussLegendre, solve
import FastGaussQuadrature  # activates Integrals' GaussLegendre node extension
using Makie: Makie, Figure, Axis, Legend, Theme, Auto, BarPlot,
             theme_latexfonts, with_theme,
             lines!, scatter!, hist!, hlines!, vlines!, band!,
             axislegend, rowsize!
using MCMCChains: MCMCChains
using Mooncake: Mooncake
using PairPlots: PairPlots, pairplot
using Printf: @printf, @sprintf
using Random: Random
using Statistics: quantile, mean, std
using Turing: Turing, @model, NUTS, MCMCThreads, sample, DynamicPPL,
              to_submodel
import FlexiChains

include("data.jl")
include("realtime.jl")
include("model.jl")
include("postprocess.jl")
include("plots.jl")
include("main.jl")

export load_linelist, build_data, bin_edges_day, which_bin, bin_labels
export prepare_model, sample_fit
export joint_model, incubation_model, transmission_delta_model,
       random_walk_rt_model
export F_offspring, F_offspring_vec
export filter_realtime, filter_by_exposure
export diagnostics, diagnostics_table, summary_table
export vector_chain, summarise, save_posterior
export plot_data, plot_rt, plot_pair, plot_posterior_predictive
export plot_delta_sense_check, plot_prior_predictives
export analyse, main

end
