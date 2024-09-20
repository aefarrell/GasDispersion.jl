using Test, Documenter, GasDispersion
using DelimitedFiles: readdlm

const GROUP = get(ENV,"GROUP","All")
const VERSION = get(ENV,"VERSION","latest")

if GROUP == "All" || GROUP == "Base"
    @testset "GasDispersion.jl tests" begin
        @test_throws MethodError scenario_builder()

        @test_throws MethodError plume()

        @test_throws MethodError puff()
    end
    # base tests
    include("base/base_types_tests.jl")

    # testing plot recipes
    include("base/recipe_tests.jl")
end

if GROUP == "All" || GROUP == "Util"
    # testing utilities
    include("utils/util_tests.jl")
end

if GROUP == "All" || GROUP == "Model"
    # testing source models
    include("source_models/jet_source_tests.jl")

    # testing dispersion models
    include("models/gaussian_plume_tests.jl")
    include("models/gaussian_puff_tests.jl")
    include("models/intpuff_tests.jl")
    include("models/simple_jet_tests.jl")
    include("models/britter_mcquaid_plume_tests.jl")
    include("models/britter_mcquaid_puff_tests.jl")
    include("models/slab_tests.jl")
end

if GROUP == "All" || GROUP == "Ext"
    #test clapeyron extension
    include("exts/clapeyron_ext_tests.jl")
end

# some doc tests don't work with julia 1.3, because of Documenter
# but I still want to run them in the tests sometimes
if GROUP == "Doc" && VERSION != "1.3"
    doctest(GasDispersion, fix=false)
end