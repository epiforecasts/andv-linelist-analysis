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

# Write the joint_model source out as a stand-alone file so the
# analysis walkthrough can pull it in via Documenter @example without
# tripping Vue's MDX parser on bare `{T}` in the Julia code.
let ll = load_linelist()
    d   = build_data(ll)
    src = @code_string joint_model(d, bin_edges_day(d.t0))
    write(joinpath(@__DIR__, "examples", "joint_model_source.jl"), src)
end

# Render Literate sources before makedocs walks the pages.
# - `docs/examples/analysis.jl` — pre-existing walkthrough; rendered with
#   `name = "analysis"` so the generated file lands as `docs/src/analysis.md`.
# - `docs/literate/*.jl`        — real-time walkthrough and any future
#   Literate-authored pages; rendered into `docs/src/`.
Literate.markdown(joinpath(@__DIR__, "examples", "analysis.jl"),
                  joinpath(@__DIR__, "src");
                  name = "analysis",
                  flavor = Literate.DocumenterFlavor(),
                  mdstrings = true, credit = false)

for path in readdir(joinpath(@__DIR__, "literate"); join = true)
    endswith(path, ".jl") || continue
    Literate.markdown(path, joinpath(@__DIR__, "src");
                      flavor = Literate.DocumenterFlavor(),
                      mdstrings = true, credit = false)
end

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
