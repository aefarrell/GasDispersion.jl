module GasDispersion

# imports
using Interpolations: Extrapolation, Line, LinearInterpolation

# exports
export Scenario, scenario_builder
export PlumeModel, Plume, plume
export GaussianPlume, SimpleJet, BritterMcQuaidPlume
export PuffModel, Puff, puff
export GaussianPuff, BritterMcQuaidPuff


# abstract types
abstract type PlumeModel end
abstract type Plume end
abstract type PuffModel end
abstract type Puff end

# helpful utilities
include("utils/scenario.jl")
include("utils/scenario_builder.jl")
include("utils/utils.jl")

"""
    plume(scenario::Scenario, model::PlumeModel)

Runs the plume dispersion model on the given scenario and returns a callable
giving the concentration of the form
    c(x, y, z[, t])

If `model` is unspecified, defaults to gaussian.

All model parameters are assumed to be in SI base units (i.e.
distances in m, velocities in m/s, mass in kg, etc.)
"""
function plume(scenario::Scenario, model::PlumeModel=GaussianPlume()) end

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
function puff(scenario::Scenario, model::PlumeModel=GaussianPuff()) end

# puff models
include("models/gaussian_puff.jl")
include("models/britter_mcquaid_puff.jl")



end
