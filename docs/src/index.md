```@meta
EditURL = "https://github.com/sbfnk/hantavirus/blob/main/README.md"
```

```@contents
Pages = ["index.md", "getting-started.md", "methods.md", "api.md"]
Depth = 2
```

```@eval
using Markdown
Markdown.parse(read(joinpath(@__DIR__, "..", "..", "README.md"), String))
```
