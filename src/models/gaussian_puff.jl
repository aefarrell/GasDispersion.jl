
# gaussian puff model
struct GaussianPuff <: PuffModel end

struct GaussianPuffSolution <: Puff
    scenario::Scenario
    model::Symbol
    downwind_dispersion::Dispersion
    crosswind_dispersion::Dispersion
    vertical_dispersion::Dispersion
end

"""
    puff(scenario::Scenario, GaussianPuff())

Generates a gaussian dispersion model on the given scenario and returns a
callable struct giving the concentration of the form
c(x, y, z, t)

"""
function puff(scenario::Scenario, model::GaussianPuff)

    required_params = [:mass_emission_rate, :release_duration, :release_height,
                       :windspeed, :pasquill_gifford]
    if all(key -> !(ismissing(getproperty(scenario,key))), required_params)
        stability = scenario.pasquill_gifford
    else
        missing_params = [ String(i) for i in filter(key -> ismissing(getproperty(scenario,key)), required_params)]
        error_string = "These parameters cannot be missing: " * join(missing_params, ", ")
        e = MissingException(error_string)
        throw(e)
    end

    # Pasquill-Gifford dispersion
    σx = crosswind_dispersion(stability; plume=false)
    σy = σx
    σz = vertical_dispersion(stability; plume=false)

    return GaussianPuffSolution(
        scenario,  #scenario::Scenario
        :gaussian, #model::Symbol
        σx, #downwind_dispersion::Dispersion
        σy, #crosswind_dispersion::Dispersion
        σz  #vertical_dispersion::Dispersion
    )

end


function (g::GaussianPuffSolution)(x,y,z,t)
    m = g.scenario.mass_emission_rate*g.scenario.release_duration
    h = g.scenario.release_height
    u = g.scenario.windspeed
    sx = g.downwind_dispersion(x)
    sy = g.crosswind_dispersion(x)
    sz = g.vertical_dispersion(x)

    c = ( m/((2*π)^(1.5)*sx*sy*sz)
        * exp(-0.5*((x-u*t)/sx)^2)
        * exp(-0.5*(y/sy)^2)
        * ( exp(-0.5*((z-h)/sz)^2) + exp(-0.5*((z+h)/sz)^2) ) )

    return c
end
