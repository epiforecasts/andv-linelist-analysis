using Test
# Use `import` (not `using`) so the package's exported `main` function is
# not brought into `Main`. Otherwise Julia's `(@main)` auto-invocation
# would run the full `analyse()` pipeline at the end of the test script.
import Hantavirus

# Exported names that are not callable functions (skip from the meta-test).
const NON_FUNCTION_EXPORTS = Symbol[]

# Functions that take no positional arguments and no keyword arguments —
# nothing to document under `# Arguments` / `# Keyword Arguments`.
function _takes_args_or_kwargs(f)
    for m in methods(f)
        # `nargs` includes the function itself, so >1 means it takes
        # at least one positional argument.
        if m.nargs > 1
            return true
        end
        # `kwarg_names` is empty for methods without kwargs.
        if !isempty(Base.kwarg_decl(m))
            return true
        end
    end
    return false
end

function _docstring_text(f)
    d = Docs.doc(f)
    return sprint(show, MIME("text/plain"), d)
end

@testset "Hantavirus.jl" begin
    @testset "exported names have docstrings with argument sections" begin
        for name in names(Hantavirus)
            name === :Hantavirus && continue
            name in NON_FUNCTION_EXPORTS && continue
            obj = getfield(Hantavirus, name)
            obj isa Function || continue
            text = _docstring_text(obj)
            # A docstring must exist — the no-docstring fallback message
            # mentions "No documentation found".
            @test !occursin("No documentation found", text)
            if _takes_args_or_kwargs(obj)
                @test occursin("Arguments", text) ||
                      occursin("Keyword Arguments", text)
            end
        end
    end
end
