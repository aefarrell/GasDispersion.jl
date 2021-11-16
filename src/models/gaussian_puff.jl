"""
    gaussian_puff_factory(scenario)

Generates a gaussian dispersion model on the given scenario and returns a
function giving the concentration of the form
c(x, y, z, t)

"""

function gaussian_puff_factory(scenario)

    required_params = [:mass_emission_rate, :release_duration, :release_height,
                       :windspeed, :pasquill_gifford]
    if all(key -> !(ismissing(getproperty(scenario,key))), required_params)

        # parameters of the jet
        m = scenario.mass_emission_rate*scenario.release_duration
        h = scenario.release_height

        # parameters of the environment
        u = scenario.windspeed
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

    function gaussian_puff(x,y,z,t)
        sx = σx(x)
        sy = σy(x)
        sz = σz(x)

        C1 = m / ( (2*π)^(1.5) * sx * sy * sz )
        C2 = exp(-0.5*((x-u*t)/sx)^2)
        C3 = exp(-0.5*(y/sy)^2)
        C4 = ( exp(-0.5*((z-h)/sz)^2) + exp(-0.5*((z+h)/sz)^2) )

        return C1*C2*C3*C4
    end

    return gaussian_puff

end
