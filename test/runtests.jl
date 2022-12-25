using GasDispersion
using Test

include("test_scenarios.jl")

@testset "GasDispersion.jl tests" begin
    @test_throws MethodError scenario_builder()

    @test_throws MethodError plume()

    @test_throws MethodError puff()
end

# testing utilities
include("util_tests.jl")

# testing plot recipes
include("recipe_tests.jl")

# testing specific models
include("gaussian_tests.jl")
include("simple_jet_tests.jl")
include("britter_tests.jl")
