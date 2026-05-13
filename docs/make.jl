using Pkg: Pkg
Pkg.instantiate()

using DocumenterVitepress
using Documenter
using Hantavirus

DocMeta.setdocmeta!(
    Hantavirus, :DocTestSetup,
    :(using Hantavirus); recursive = true
)

# Stage repo-root figures as vitepress public assets so the README's
# `figures/...` image links resolve on the docs site (vitepress otherwise
# treats them as rollup module imports and fails the build).
cp(joinpath(@__DIR__, "..", "figures"),
   joinpath(@__DIR__, "src", "public", "figures"); force = true)

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

deploydocs(;
    repo = "github.com/sbfnk/hantavirus",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
