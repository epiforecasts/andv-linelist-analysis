```@meta
EditURL = "https://github.com/sbfnk/hantavirus/blob/main/README.md"
```

```@eval
using Markdown
raw = read(joinpath(@__DIR__, "..", "..", "README.md"), String)
# Adapt repo-root-relative links in the README for the docs site:
# - figures/* are staged into docs/src/public/figures/ by make.jl, so reference
#   them via the absolute public path /figures/...
# - METHODS.md has its own docs page, so link there instead of the markdown file
# - LICENSE has no docs page; drop the link (text remains).
rewritten = replace(raw,
    r"\]\(figures/" => "](/figures/",
    "](METHODS.md)" => "](methods)",
    "[LICENSE](LICENSE)" => "LICENSE",
)
Markdown.parse(rewritten)
```
