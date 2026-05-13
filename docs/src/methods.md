```@meta
EditURL = "https://github.com/sbfnk/hantavirus/blob/main/METHODS.md"
```

```@eval
using Markdown
raw = read(joinpath(@__DIR__, "..", "..", "METHODS.md"), String)
# METHODS.md links back to README.md (a repo-relative path); on the docs site
# that's the home page at /.
rewritten = replace(raw, "[README](README.md)" => "[README](/)")
Markdown.parse(rewritten)
```
