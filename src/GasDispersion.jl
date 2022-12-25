module GasDispersion

# imports
using Markdown
using Interpolations: Extrapolation, Line, LinearInterpolation
using SpecialFunctions: erf
using RecipesBase: @recipe

# source models
export Atmosphere, Ambient
export Release, Scenario, scenario_builder
export SourceModel, JetSource

# plume models
export PlumeModel, Plume, plume
export GaussianPlume, SimpleJet, BritterMcQuaidPlume

# puff models
export PuffModel, Puff, puff
export GaussianPuff, IntPuff, BritterMcQuaidPuff


# abstract types
abstract type Atmosphere end
abstract type SourceModel end
abstract type PlumeModel end
abstract type Plume end
abstract type PuffModel end
abstract type Puff end

# helpful utilities
include("utils/scenario.jl")
include("utils/atmosphere.jl")
include("utils/utils.jl")

# plot recipes
include("recipes/plume_recipes.jl")


"""
    plume(scenario::Scenario, model::PlumeModel)

Runs the plume dispersion model on the given scenario and returns a callable
giving the concentration of the form
    c(x, y, z[, t])

If `model` is unspecified, defaults to gaussian.

All model parameters are assumed to be in SI base units (i.e.
distances in m, velocities in m/s, mass in kg, etc.)
"""
function plume end

#default behaviour
plume(s) = plume(s, GaussianPlume())

# plume models
include("models/gaussian_plume.jl")
include("models/simple_jet.jl")
include("models/britter_mcquaid_plume.jl")


"""
    puff(scenario::Scenario, model::PuffModel)

Runs the puff dispersion model on the given scenario and returns a callable
giving the concentration of the form
    c(x, y, z, t)

If `model` is unspecified, defaults to gaussian.

All model parameters are assumed to be in SI base units (i.e.
distances in m, velocities in m/s, mass in kg, etc.)
"""
function puff end

#default behaviour
puff(s) = puff(s, GaussianPuff())

# puff models
include("models/gaussian_puff.jl")
include("models/intpuff.jl")
include("models/britter_mcquaid_puff.jl")


"""
    scenario_builder(source::SourceModel, atmosphere::Atmosphere)

Builds a scenario given a source model and an atmosphere. If no atmosphere is
given defaults to ambient conditions.

"""
function scenario_builder end

#default behaviour
scenario_builder(m) = scenario_builder(m, Ambient())

# source models
include("source_models/jet_source.jl")


end
