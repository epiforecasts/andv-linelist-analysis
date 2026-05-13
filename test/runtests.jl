using Test
using Aqua
using Hantavirus

@testset "Aqua.jl meta-tests" begin
    @testset "Unbound args" begin
        Aqua.test_unbound_args(Hantavirus)
    end

    @testset "Undefined exports" begin
        Aqua.test_undefined_exports(Hantavirus)
    end

    @testset "Project extras" begin
        Aqua.test_project_extras(Hantavirus)
    end

    @testset "Stale deps" begin
        # CairoMakie is listed as a dep so users get a working Makie
        # backend when they call plotting helpers; the package itself
        # never does `using CairoMakie`, so Aqua flags it as stale.
        Aqua.test_stale_deps(Hantavirus; ignore = [:CairoMakie])
    end

    @testset "Deps compat" begin
        Aqua.test_deps_compat(Hantavirus)
    end

    @testset "Undocumented names" begin
        Aqua.test_undocumented_names(Hantavirus)
    end

    @testset "Piracies" begin
        Aqua.test_piracies(Hantavirus)
    end

    @testset "Ambiguities" begin
        Aqua.test_ambiguities(Hantavirus)
    end
end
