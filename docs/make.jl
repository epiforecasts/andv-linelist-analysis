using Pkg: Pkg
Pkg.instantiate()

using Documenter
using Hantavirus

DocMeta.setdocmeta!(
    Hantavirus, :DocTestSetup,
    :(using Hantavirus); recursive = true
)

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
        "Getting Started" => "getting-started.md",
        "Methods" => "methods.md",
        "API Reference" => "api.md",
    ],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://sbfnk.github.io/hantavirus",
        repolink = "https://github.com/sbfnk/hantavirus",
    ),
)

deploydocs(;
    repo = "github.com/sbfnk/hantavirus",
    devbranch = "main",
    push_preview = true,
)
