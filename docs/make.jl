using Pkg: Pkg
Pkg.instantiate()

using DocumenterVitepress
using Documenter
using Literate
using Hantavirus

DocMeta.setdocmeta!(
    Hantavirus, :DocTestSetup,
    :(using Hantavirus); recursive = true
)

# Render Literate sources into docs/src before makedocs walks the pages.
# `execute = false` defers code execution to Documenter's @example blocks,
# so the rendered markdown lives in the same session as other pages.
const LITERATE_DIR  = joinpath(@__DIR__, "literate")
const GENERATED_DIR = joinpath(@__DIR__, "src")
for path in readdir(LITERATE_DIR; join = true)
    endswith(path, ".jl") || continue
    Literate.markdown(path, GENERATED_DIR;
                      documenter = true, execute = false)
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
        "Methods" => "methods.md",
        "Real-time monitoring" => "realtime.md",
        "API Reference" => "api.md",
    ],
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/sbfnk/hantavirus",
        devbranch = "main",
        devurl = "dev",
    ),
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/sbfnk/hantavirus",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
