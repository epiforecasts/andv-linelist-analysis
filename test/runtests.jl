using Test
using Aqua
using TransmissionLinelist

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
