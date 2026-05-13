```@meta
EditURL = "https://github.com/sbfnk/hantavirus/blob/main/README.md"
```

```@eval
using Markdown
raw = read(joinpath(@__DIR__, "..", "..", "README.md"), String)
# Repo-root-relative image paths in the README (e.g. `figures/Rt.png`) need to
# become vitepress public-asset paths (`/figures/Rt.png`) for the docs build.
# The corresponding files are staged into docs/src/public/figures/ by make.jl.
rewritten = replace(raw, r"\]\(figures/" => "](/figures/")
Markdown.parse(rewritten)
```
