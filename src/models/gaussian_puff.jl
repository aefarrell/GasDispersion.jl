
# gaussian puff model
struct GaussianPuff <: PuffModel end

struct GaussianPuffSolution <: Puff
    scenario::Scenario
    model::Symbol
    mass::Number
    height::Number
    windspeed::Number
    downwind_dispersion::Dispersion
    crosswind_dispersion::Dispersion
    vertical_dispersion::Dispersion
end

@doc doc"""
    puff(scenario::Scenario, GaussianPuff())

Generates a gaussian dispersion model on the given scenario and returns a
callable struct giving the concentration of the form
c(x, y, z, t)

```math
c\left(x,y,z,t\right) = { m \over { (2 \pi)^{3/2} \sigma_x \sigma_y \sigma_z } }
\exp \left( -\frac{1}{2} \left( {x - ut} \over \sigma_x \right)^2 \right)
\exp \left( -\frac{1}{2} \left( {y} \over \sigma_y \right)^2 \right)
\times \left[ \exp \left( -\frac{1}{2} \left( {z - h} \over \sigma_z \right)^2 \right)
+ \exp \left( -\frac{1}{2} \left( {z + h} \over \sigma_z \right)^2 \right)\right]
```

"""
function puff(scenario::Scenario, model::GaussianPuff)

    stability = scenario.atmosphere.stability
    m = scenario.release.mass_rate*scenario.release.duration
    h = scenario.release.height
    u = scenario.atmosphere.windspeed

    # Pasquill-Gifford dispersion
    σx = crosswind_dispersion(stability; plume=false)
    σy = σx
    σz = vertical_dispersion(stability; plume=false)

    return GaussianPuffSolution(
        scenario,  #scenario::Scenario
        :gaussian, #model::Symbol
        m,  #mass
        h,  #release height
        u,  #windspeed
        σx, #downwind_dispersion::Dispersion
        σy, #crosswind_dispersion::Dispersion
        σz  #vertical_dispersion::Dispersion
    )

end


function (g::GaussianPuffSolution)(x,y,z,t)
    m = g.mass
    h = g.height
    u = g.windspeed
    sx = g.downwind_dispersion(x)
    sy = g.crosswind_dispersion(x)
    sz = g.vertical_dispersion(x)

    c = ( m/((2*π)^(1.5)*sx*sy*sz)
        * exp(-0.5*((x-u*t)/sx)^2)
        * exp(-0.5*(y/sy)^2)
        * ( exp(-0.5*((z-h)/sz)^2) + exp(-0.5*((z+h)/sz)^2) ) )

    return c
end
