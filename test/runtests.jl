using Test, Documenter, GasDispersion
using DelimitedFiles: readdlm

const GROUP = get(ENV,"GROUP","All")
const VERSION = get(ENV,"VERSION","latest")

if GROUP ∈ ["All", "Base"]
    @testset "GasDispersion.jl tests" begin
        @test_throws MethodError scenario_builder()

        @test_throws MethodError plume()

        @test_throws MethodError puff()
    end
    # base tests
    include("base/base_types_tests.jl")

    # testing plot recipes
    include("base/recipe_tests.jl")

    # testing utilities
    include("utils/util_tests.jl")
    include("depreciation/depreciation_tests.jl")

    # testing source models
    include("source_models/jet_source_tests.jl")

    # testing dispersion models
    include("models/gaussian_plume_tests.jl")
    include("models/gaussian_ml_tests.jl")
    include("models/gaussian_puff_tests.jl")
    include("models/intpuff_tests.jl")
    include("models/palazzi_tests.jl")
    include("models/simple_jet_tests.jl")
    include("models/britter_mcquaid_plume_tests.jl")
    include("models/britter_mcquaid_puff_tests.jl")
    include("models/slab_tests.jl")
end

if GROUP ∈ ["All", "Ext"] && VERSION ∈ ["1","latest","nightly"]
    # this feels kind of janky to me, but it stops versions <1.9
    # from trying to run the Clapeyron extension
    import Pkg; Pkg.add("Clapeyron")
    using Clapeyron
    include("exts/clapeyron_ext_tests.jl")
end

# some doc tests don't work with julia 1.3, because of Documenter
# but I still want to run them in the tests sometimes
if GROUP ∈ ["All", "Doc"] && VERSION != "1.3"
    doctest(GasDispersion, fix=false)
end