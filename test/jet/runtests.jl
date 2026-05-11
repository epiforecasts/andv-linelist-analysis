using Test
using JET
using Hantavirus

# JET reports false positives inside Turing @model blocks: the `~` macro
# rewrites variable assignments into varinfo lookups that JET cannot see
# through. All six current reports are of this kind (undefined locals in
# joint_model). Marked @test_broken until JET or DynamicPPL gains support
# for analysing @model-generated code.
@test_broken JET.test_package(
    Hantavirus; target_modules = (Hantavirus,)
)
