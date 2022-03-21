module GasDispersion

export Scenario, Plume, Puff, plume, puff, scenario_builder

using Interpolations

# helpful utilities
include("utils/scenario.jl")
include("utils/scenario_builder.jl")
include("utils/utils.jl")

# plume models
include("models/gaussian_plume.jl")
include("models/britter_mcquaid_plume.jl")

# puff models
include("models/gaussian_puff.jl")
include("models/britter_mcquaid_puff.jl")

"""
    plume(scenario::Scenario; model=:gaussian, kwargs...)

Runs the plume dispersion model on the given scenario and returns a function
giving the concentration of the form
    c(x, y, z[, t])

If `model` is unspecified, defaults to gaussian, `kwargs` are passed to the
plume model. All model parameters are assumed to be in SI base units (i.e.
distances in m, velocities in m/s, mass in kg, etc.)
"""
function plume(scenario::Scenario; model=:gaussian, kwargs...)
    if model==:gaussian
        return gaussian_plume_factory(scenario; kwargs...)
    elseif model==:brittermcquaid
        return britter_plume_factory(scenario; kwargs...)
    else
        error_string = "plume dispersion model $model is not currently implemented"
        error(error_string)
    end
end

"""
    puff(scenario::Scenario; model=:gaussian, kwargs...)

Runs the puff dispersion model on the given scenario and returns a function
giving the concentration of the form
    c(x, y, z, t)

If `model` is unspecified, defaults to gaussian, `kwargs` are passed to the
puff model. All model parameters are assumed to be in SI base units (i.e.
distances in m, velocities in m/s, mass in kg, etc.)
"""
function puff(scenario::Scenario; model=:gaussian, kwargs...)
    if model==:gaussian
        return gaussian_puff_factory(scenario; kwargs...)
    # elseif model=="britter-mcquaid"
    #     return britter_plume_factory(scenario; kwargs...)
    else
        error_string = "puff dispersion model $model is not currently implemented"
        error(error_string)
    end
end

end
