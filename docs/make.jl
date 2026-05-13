using Pkg: Pkg
Pkg.instantiate()

using DocumenterVitepress
using Documenter
using Hantavirus

DocMeta.setdocmeta!(
    Hantavirus, :DocTestSetup,
    :(using Hantavirus); recursive = true
)

makedocs(;
    sitename = "Hantavirus.jl",
    authors = "Sebastian Funk, Sam Abbott, and contributors",
    clean = true,
    doctest = false,
    linkcheck = true,
    warnonly = [:docs_block, :missing_docs, :autodocs_block, :linkcheck],
    modules = [Hantavirus],
    pages = [
        "Home" => "index.md",
        "Methods" => "methods.md",
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
