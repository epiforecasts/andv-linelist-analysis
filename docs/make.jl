using Pkg: Pkg
Pkg.instantiate()

using DocumenterVitepress
using Documenter
using Literate
using CodeTracking
using CairoMakie
using TransmissionLinelist

# Retina-quality figures in the rendered docs.
CairoMakie.activate!(; px_per_unit = 2.0)

DocMeta.setdocmeta!(
    TransmissionLinelist, :DocTestSetup,
    :(using TransmissionLinelist); recursive = true
)

# Write the joint model `@model` definition source out as a stand-alone
# fenced code block so the Literate page can pull it in via Documenter
# @example / explicit include without going through Vue's MDX parser,
# which trips on bare `{T}` in the Julia code.
let ll = load_linelist()
    d = build_data(ll)
    src = @code_string TransmissionLinelist.joint_model_def(d, bin_edges_day(d.t0))
    write(joinpath(@__DIR__, "examples", "joint_model_source.jl"), src)
end

# Render the executable walkthroughs through Literate so figures and tables
# are generated at build time from a single source script each.
const LITERATE_OUT = joinpath(@__DIR__, "src")

Literate.markdown(joinpath(@__DIR__, "examples", "analysis.jl"),
    LITERATE_OUT;
    name = "analysis",
    flavor = Literate.DocumenterFlavor(),
    mdstrings = true, credit = false)

Literate.markdown(joinpath(@__DIR__, "examples", "realtime.jl"),
    LITERATE_OUT;
    name = "realtime",
    flavor = Literate.DocumenterFlavor(),
    mdstrings = true, credit = false)

makedocs(;
    # Short label for the rendered Vitepress sidebar header and browser tab.
    # The full descriptive title lives as the h1 of `docs/src/index.md`, so
    # the rendered landing page still leads with the long form.
    sitename = "ANDV linelist analysis",
    authors = "Sebastian Funk, Sam Abbott, and contributors",
    clean = true,
    doctest = false,
    linkcheck = true,
    warnonly = [:docs_block, :missing_docs, :autodocs_block, :linkcheck],
    modules = [TransmissionLinelist],
    pages = [
        "Home" => "index.md",
        "Model" => "model.md",
        "Limitations" => "limitations.md",
        "Analysis walkthrough" => "analysis.md",
        "Real-time monitoring" => "realtime.md",
        "API Reference" => "api.md"
    ],
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/epiforecasts/andv-linelist-analysis",
        devbranch = "main",
        devurl = "dev"
    )
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/epiforecasts/andv-linelist-analysis",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true
)
