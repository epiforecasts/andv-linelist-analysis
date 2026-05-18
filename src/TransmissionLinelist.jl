module TransmissionLinelist

using ADTypes: AutoMooncake
using AlgebraOfGraphics: AlgebraOfGraphics, data, mapping, visual, draw!,
                         nonnumeric
using ArgParse: ArgParseSettings, @add_arg_table!, parse_args
using Chain: @chain
using CSV: CSV
using DataFrames: DataFrame, nrow, eachrow, passmissing, rename!, sort!
using DataFramesMeta: @select, @transform, @subset, @combine, @by, @rtransform,
                      @rsubset, @orderby, @rename
using Dates: Dates, Date, Day, dayofweek
using Distributions: Distributions, Normal, LogNormal, truncated,
                     NegativeBinomial, Uniform, Gamma, Poisson,
                     ContinuousUnivariateDistribution,
                     logpdf, logcdf, cdf, ccdf, pdf, partype
using DocStringExtensions: TYPEDSIGNATURES
using Integrals: IntegralProblem, GaussLegendre, solve
import FastGaussQuadrature  # activates Integrals' GaussLegendre node extension
using Mooncake: Mooncake
using PairPlots: PairPlots, pairplot
using Makie: Makie, Figure, Axis, Legend, Theme, Auto, BarPlot, Hist,
             Lines, Scatter, MarkerElement, PolyElement, GridLayout,
             theme_latexfonts, with_theme,
             lines!, scatter!, hist!, hlines!, vlines!, band!,
             barplot!, rangebars!,
             axislegend, rowsize!, colsize!
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

export load_linelist, build_data, prepare_rt_edges, log_R_at, bin_labels
export joint_model, delays_only_model
export incubation_model, transmission_delta_model, random_walk_rt_model,
       nb_dispersion_model, truncation_model, case_model,
       latent_times_model, posterior_predictive
export ConvolvedDelays
export filter_realtime, filter_by_exposure
export predict_controlled_outbreak, predict_natural_chain_outbreak,
       realised_future_count
export diagnostics, diagnostics_table, summary_table
export vector_chain, summarise, save_posterior, rt_band, summarise_predictive
export plot_data, plot_rt, plot_pair, plot_predictive_distributions
export plot_delta_sense_check, plot_inc_sense_check, plot_z_ppc
export plot_prior_predictives
export z_ppc_summary
export analyse, sample_fit, main

end
