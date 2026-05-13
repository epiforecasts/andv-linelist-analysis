using Pkg: Pkg
Pkg.instantiate()

using DocumenterVitepress
using Documenter
using TransmissionLinelist

DocMeta.setdocmeta!(
    TransmissionLinelist, :DocTestSetup,
    :(using TransmissionLinelist); recursive = true
)

makedocs(;
    sitename = "TransmissionLinelist.jl",
    authors = "Sebastian Funk, Sam Abbott, and contributors",
    clean = true,
    doctest = false,
    linkcheck = true,
    warnonly = [:docs_block, :missing_docs, :autodocs_block, :linkcheck],
    modules = [TransmissionLinelist],
    pages = [
        "Home" => "index.md",
        "Methods" => "methods.md",
        "API Reference" => "api.md",
    ],
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/sbfnk/andv-linelist-analysis",
        devbranch = "main",
        devurl = "dev",
    ),
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/sbfnk/andv-linelist-analysis",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
