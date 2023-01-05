module GasDispersion

# imports
using Markdown
using Interpolations: Extrapolation, Line, LinearInterpolation
using SpecialFunctions: erf
using RecipesBase

# source models
export Atmosphere, DryAir
export StabilityClass, ClassA, ClassB, ClassC, ClassD, ClassE, ClassF
export Substance, Release, Scenario, scenario_builder
export SourceModel, JetSource

# plume models
export PlumeModel, Plume, plume
export GaussianPlume, SimpleJet, BritterMcQuaidPlume

# puff models
export PuffModel, Puff, puff
export GaussianPuff, IntPuff, BritterMcQuaidPuff


# abstract types
abstract type Atmosphere end
abstract type StabilityClass end
abstract type SourceModel end
abstract type PlumeModel end
abstract type Plume end
abstract type PuffModel end
abstract type Puff end

# basic type definitions and such
include("base/base_types.jl")
include("base/plot_recipes.jl")

# helpful utilities
include("utils/utils.jl")


"""
    plume(scenario::Scenario, model::PlumeModel)

Runs the plume dispersion model on the given scenario and returns a callable
giving the concentration of the form
    c(x, y, z[, t])

The concentration is in kg/m³, if `model` is unspecified, defaults to gaussian.

All model parameters are assumed to be in SI base units (i.e.
distances in m, velocities in m/s, mass in kg, etc.)
"""
function plume end

#default behaviour
plume(s; kwargs...) = plume(s, GaussianPlume; kwargs...)

# plume models
include("models/gaussian_plume.jl")
include("models/simple_jet.jl")
include("models/britter_mcquaid_plume.jl")


"""
    puff(scenario::Scenario, model::PuffModel)

Runs the puff dispersion model on the given scenario and returns a callable
giving the concentration of the form
    c(x, y, z, t)

The concentration is in kg/m³, if `model` is unspecified, defaults to gaussian.

All model parameters are assumed to be in SI base units (i.e.
distances in m, velocities in m/s, mass in kg, etc.)
"""
function puff end

# default behaviour
puff(s; kwargs...) = puff(s, GaussianPuff; kwargs...)

# puff models
include("models/gaussian_puff.jl")
include("models/intpuff.jl")
include("models/britter_mcquaid_puff.jl")


"""
    scenario_builder(substance::Substance, source::SourceModel, atmosphere::Atmosphere)

Builds a scenario given a substance, source model and an atmosphere.
If no atmosphere is given defaults to ambient conditions.

"""
function scenario_builder end

# default behaviour
scenario_builder(s,m; kwargs...) = scenario_builder(s,m,DryAir(); kwargs...)

# source models
include("source_models/jet_source.jl")


end
