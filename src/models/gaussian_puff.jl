"""
    gaussian_puff_factory(scenario)

Generates a gaussian dispersion model on the given scenario and returns a
function giving the concentration of the form
c(x, y, z, t)

"""

function gaussian_puff_factory(scenario)
    # parameters of the jet
    m = scenario.mass_emission_rate
    Dⱼ = scenario.jet_diameter
    uⱼ = scenario.jet_velocity
    h = scenario.release_height

    # parameters of the environment
    u = scenario.windspeed
    stability = scenario.pasquill_gifford

    # Pasquill-Gifford dispersion
    σx = crosswind_dispersion(stability; plume=false)
    σy = crosswind_dispersion(stability; plume=false)
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
