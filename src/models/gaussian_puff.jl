
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
callable giving the concentration of the form `c(x, y, z, t)`

Gaussian puff model is per *Guidelines for Consequence Analysis of Chemical
Release*, CCPS, New York (1999)

```math
c\left(x,y,z,t\right) = { {\dot{m} \Delta t} \over n }
{ { \exp \left( -\frac{1}{2} \left( {x - u t } \over \sigma_x \right)^2 \right) } \over { \sqrt{2\pi} \sigma_x } }
{ { \exp \left( -\frac{1}{2} \left( {y} \over \sigma_y \right)^2 \right) } \over { \sqrt{2\pi} \sigma_y } }\\
\times { { \exp \left( -\frac{1}{2} \left( {z - h} \over \sigma_z \right)^2 \right)
+ \exp \left( -\frac{1}{2} \left( {z + h} \over \sigma_z \right)^2 \right) } \over { \sqrt{2\pi} \sigma_z } }
```

"""
function puff(scenario::Scenario, model::GaussianPuff)

    stability = scenario.atmosphere.stability
    G = scenario.release.mass_rate*scenario.release.duration
    h = scenario.release.height
    u = scenario.atmosphere.windspeed

    # Pasquill-Gifford dispersion
    σx = crosswind_dispersion(stability; plume=false)
    σy = σx
    σz = vertical_dispersion(stability; plume=false)

    return GaussianPuffSolution(
        scenario,  #scenario::Scenario
        :gaussian, #model::Symbol
        G,  #mass
        h,  #release height
        u,  #windspeed
        σx, #downwind_dispersion::Dispersion
        σy, #crosswind_dispersion::Dispersion
        σz  #vertical_dispersion::Dispersion
    )

end


function (g::GaussianPuffSolution)(x,y,z,t)

    # domain check
    if (x<0)||(z<0)||(t<0)
        return 0.0
    end

    G = g.mass
    h = g.height
    u = g.windspeed
    σx = g.downwind_dispersion(x)
    σy = g.crosswind_dispersion(x)
    σz = g.vertical_dispersion(x)

    c = ( G/((2*π)^(1.5)*σx*σy*σz)
        * exp(-0.5*((x-u*t)/σx)^2)
        * exp(-0.5*(y/σy)^2)
        * ( exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2) ) )

    return c
end
