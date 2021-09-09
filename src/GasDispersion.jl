module GasDispersion

export Scenario, plume

# helpful utilities
include("utils/scenario_builder.jl")
include("utils/utils.jl")

# plume models
include("models/gaussian_plume.jl")

# puff models


struct Scenario
    mass_emission_rate::Number
    jet_diameter::Number
    jet_velocity::Number
    jet_density::Number
    release_pressure::Number
    release_temperature::Number
    release_height::Number
    windspeed::Number
    ambient_density::Number
    ambient_pressure::Number
    ambient_temperature::Number
    pasquill_gifford::String
end

"""
    plume(scenario::Scenario; model="gaussian", kwargs...)

Runs the plume dispersion model on the given scenario and returns a function
giving the concentration of the form
    c(x, y, z[, t])

If `model` is unspecified, defaults to gaussian, `kwargs` are passed to the
plume model. All model parameters are assumed to be in SI base units (i.e.
distances in m, velocities in m/s, mass in kg, etc.)
"""
function plume(scenario::Scenario; model::String="gaussian", kwargs...)
    if model=="gaussian"
        conc_fun= gaussian_plume_factory(scenario; kwargs...)
    else
        error_string = string("plume dispersion model ''",model,"'' is not currently implemented")
        error(error_string)
    end

    return conc_fun
end

end
