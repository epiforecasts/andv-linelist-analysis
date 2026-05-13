```@meta
EditURL = "https://github.com/sbfnk/hantavirus/blob/main/README.md"
```

```@eval
using Markdown
raw = read(joinpath(@__DIR__, "..", "..", "README.md"), String)
# Rewrite repo-root-relative image paths to absolute raw URLs so the rendered
# docs site (built from docs/) can resolve them.
rewritten = replace(raw,
    r"\]\(figures/" => "](https://raw.githubusercontent.com/sbfnk/hantavirus/main/figures/")
Markdown.parse(rewritten)
```
