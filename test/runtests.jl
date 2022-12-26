using GasDispersion
using Test

@testset "GasDispersion.jl tests" begin
    @test_throws MethodError scenario_builder()

    @test_throws MethodError plume()

    @test_throws MethodError puff()
end

# testing scenarios
include("test_scenarios.jl")

# testing utilities
include("util_tests.jl")

# testing specific models
include("gaussian_tests.jl")
include("simple_jet_tests.jl")
include("britter_tests.jl")

# testing plot recipes
# this is comedically slow
include("recipe_tests.jl")