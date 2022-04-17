using GasDispersion
using Test

include("test_scenarios.jl")

@testset "GasDispersion.jl tests" begin
    @test_throws MethodError Scenario()

    @test_throws MethodError plume()
    @test_throws ErrorException plume(test_scenario, model=:error)
    @test_throws ErrorException plume(bad_class)

    @test_throws MethodError puff()
    @test_throws ErrorException puff(test_scenario, model=:error)
    @test_throws ErrorException puff(bad_class)

end

# testing utilities
include("util_tests.jl")

# testing specific models
include("gaussian_tests.jl")
include("simple_jet_tests.jl")
include("britter_tests.jl")
