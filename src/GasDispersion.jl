__precompile__()

module GasDispersion

# imports
using Markdown
using DataInterpolations: LinearInterpolation, AkimaInterpolation
using SpecialFunctions: erf
using RecipesBase

# releases
export Atmosphere, SimpleAtmosphere
export StabilityClass, ClassA, ClassB, ClassC, ClassD, ClassE, ClassF
export AbstractSubstance, Substance
export Release, HorizontalJet, VerticalJet
export Scenario, scenario_builder

# source models
export SourceModel, JetSource

# plume models
export PlumeModel, Plume, plume
export GaussianPlume, GaussianMixingLayer, SimpleJet, BritterMcQuaidPlume

# puff models
export PuffModel, Puff, puff
export GaussianPuff, Palazzi, IntPuff, BritterMcQuaidPuff, SLAB

# equation sets
export EquationSet, BasicEquationSet, DefaultSet, DefaultPuffSet
export CCPSRural, CCPSUrban, TNOPlume, Turner, ISC3Rural, ISC3Urban
export CCPSPuffRural, CCPSPuffUrban, TNOPuff

# windspeed correlations
export DefaultWind,IrwinRural, IrwinUrban, ISC3UrbanWind
export BusingerWind, TNOWind

# dispersion correlations
export Defaultσy, Defaultσz
export BriggsRuralσy, BriggsUrbanσy, ISC3Ruralσy, TNOPlumeσy, Turnerσy
export BriggsRuralσz, BriggsUrbanσz, ISC3Ruralσz, TNOPlumeσz, Turnerσz
export CCPSPuffσx, CCPSPuffσy, CCPSPuffσz
export TNOPuffσx, TNOPuffσy, TNOPuffσz

# abstract types
abstract type AbstractSubstance end
abstract type Atmosphere end
abstract type StabilityClass end
abstract type Release end
abstract type SourceModel end
abstract type PlumeModel end
abstract type Plume end
abstract type PuffModel end
abstract type Puff end
abstract type EquationSet end
abstract type PowerLawWind end
abstract type MoninObukhovWind end
abstract type DispersionFunction end

const R_GAS_CONST = 8.31446261815324 #Universal Gas Constant, J/mol/K

# basic type definitions and such
include("base/base_types.jl")
include("base/plot_recipes.jl")

# helpful utilities
include("utils/utils.jl")


"""
    plume(scenario::Scenario, model::PlumeModel[, equationset::EquationSet])

Runs the plume dispersion model on the given scenario and returns the solution
which is callable to give the concentration
    c(x, y, z[, t])

The concentration is in vol fraction, if `model` is unspecified, defaults to a simple
gaussian plume model.

`equationset`s are used to specify that an alternative set of correlations 
should be used for model parameters, if alternatives exist.

All model parameters are assumed to be in SI base units (i.e.
distances in m, velocities in m/s, mass in kg, etc.)
"""
function plume end

#default behaviour
plume(s; kwargs...) = plume(s, GaussianPlume(); kwargs...)

# plume models
include("models/gaussian_plume.jl")
include("models/gaussian_mixing_layer.jl")
include("models/simple_jet.jl")
include("models/britter_mcquaid_plume.jl")


"""
    puff(scenario::Scenario, model::PuffModel[, equationset::EquationSet])

Runs the puff dispersion model on the given scenario and returns the solution
which is callable to give the concentration
    c(x, y, z, t)

The concentration is in vol fraction, if `model` is unspecified, defaults to a
simple gaussian puff.

`equationset`s are used to specify that an alternative set of correlations 
should be used for model parameters, if alternatives exist.

All model parameters are assumed to be in SI base units (i.e.
distances in m, velocities in m/s, mass in kg, etc.)
"""
function puff end

# default behaviour
puff(s; kwargs...) = puff(s, GaussianPuff(); kwargs...)

# puff models
include("models/gaussian_puff.jl")
include("models/palazzi_puff.jl")
include("models/intpuff.jl")
include("models/britter_mcquaid_puff.jl")
include("models/slab_puff.jl")


"""
    scenario_builder(substance::Substance, source::SourceModel, atmosphere::Atmosphere)

Builds a scenario given a substance, source model and an atmosphere.
If no atmosphere is given defaults to dry air at ambient conditions and class
F stability.

"""
function scenario_builder end

# default behaviour
scenario_builder(s,m; kwargs...) = scenario_builder(s,m,SimpleAtmosphere(); kwargs...)

# source models
include("source_models/jet_source.jl")


end
