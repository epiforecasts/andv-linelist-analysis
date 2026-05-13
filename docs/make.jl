using Pkg: Pkg
Pkg.instantiate()

using DocumenterVitepress
using Documenter
using Literate
using CodeTracking
using CairoMakie
using Hantavirus

# Retina-quality figures in the rendered docs.
CairoMakie.activate!(; px_per_unit = 2.0)

DocMeta.setdocmeta!(
    Hantavirus, :DocTestSetup,
    :(using Hantavirus); recursive = true
)

let ll = load_linelist()
    d   = build_data(ll)
    src = @code_string Hantavirus.joint_model_def(d, bin_edges_day(d.t0))
    write(joinpath(@__DIR__, "examples", "joint_model_source.jl"), src)
end

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
    sitename = "Hantavirus.jl",
    authors = "Sebastian Funk, and contributors",
    clean = true,
    doctest = false,
    linkcheck = true,
    warnonly = [:docs_block, :missing_docs, :autodocs_block, :linkcheck],
    modules = [Hantavirus],
    pages = [
        "Home" => "index.md",
        "Model" => "model.md",
        "Limitations" => "limitations.md",
        "Analysis walkthrough" => "analysis.md",
        "Real-time monitoring" => "realtime.md",
        "API Reference" => "api.md",
    ],
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/sbfnk/hantavirus",
        devbranch = "main",
        devurl = "dev",
    ),
)

deploydocs(;
    repo = "github.com/sbfnk/hantavirus",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
