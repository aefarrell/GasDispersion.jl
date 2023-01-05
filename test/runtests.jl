using GasDispersion
using Test

@testset "GasDispersion.jl tests" begin
    @test_throws MethodError scenario_builder()

    @test_throws MethodError plume()

    @test_throws MethodError puff()
end

# base tests
include("base/base_types_tests.jl")

# testing utilities
include("utils/util_tests.jl")

# # testing source models
# include("source_models/jet_source_tests.jl")
#
# # testing dispersion models
# include("models/gaussian_tests.jl")
# include("models/simple_jet_tests.jl")
# include("models/britter_tests.jl")
#
# # testing plot recipes
# # this is comedically slow
# include("base/recipe_tests.jl")
