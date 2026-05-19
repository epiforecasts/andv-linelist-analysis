using Test
using Aqua
# Use `import` so the package's exported `main` does not shadow
# downstream uses inside the testset bodies.
import TransmissionLinelist
using TransmissionLinelist: ConvolvedDelays
using DifferentiationInterface: AutoMooncake, AutoForwardDiff,
                                value_and_gradient
using Mooncake: Mooncake

include("test_realtime.jl")
include("test_recovery.jl")
include("test_submodel_recovery.jl")
include("test_jet.jl")

@testset "Aqua.jl meta-tests" begin
    @testset "Unbound args" begin
        Aqua.test_unbound_args(TransmissionLinelist)
    end

    @testset "Undefined exports" begin
        Aqua.test_undefined_exports(TransmissionLinelist)
    end

    @testset "Project extras" begin
        Aqua.test_project_extras(TransmissionLinelist)
    end

    @testset "Stale deps" begin
        Aqua.test_stale_deps(TransmissionLinelist)
    end

    @testset "Deps compat" begin
        Aqua.test_deps_compat(TransmissionLinelist)
    end

    @testset "Undocumented names" begin
        Aqua.test_undocumented_names(TransmissionLinelist)
    end

    @testset "Piracies" begin
        Aqua.test_piracies(TransmissionLinelist)
    end

    @testset "Ambiguities" begin
        Aqua.test_ambiguities(TransmissionLinelist)
    end
end

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

@testset "TransmissionLinelist.jl" begin
    @testset "exported names have docstrings with argument sections" begin
        for name in names(TransmissionLinelist)
            name === :TransmissionLinelist && continue
            name in NON_FUNCTION_EXPORTS && continue
            obj = getfield(TransmissionLinelist, name)
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
