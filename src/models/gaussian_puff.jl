
# gaussian puff model
struct GaussianPuff <: PuffModel end

struct GaussianPuffSolution{S<:StabilityClass} <: Puff
    scenario::Scenario
    model::Symbol
    mass::Number
    height::Number
    windspeed::Number
    stability::Type{S}
end

@doc doc"""
    puff(scenario::Scenario, GaussianPuff)

Returns the solution to a Gaussian puff dispersion model for the given scenario.

The Gaussian puff model is per *Guidelines for Consequence Analysis of Chemical
Release*, CCPS, New York (1999)

```math
c\left(x,y,z,t\right) = \dot{m} \Delta t
{ { \exp \left( -\frac{1}{2} \left( {x - u t } \over \sigma_x \right)^2 \right) } \over { \sqrt{2\pi} \sigma_x } }
{ { \exp \left( -\frac{1}{2} \left( {y} \over \sigma_y \right)^2 \right) } \over { \sqrt{2\pi} \sigma_y } }\\
\times { { \exp \left( -\frac{1}{2} \left( {z - h} \over \sigma_z \right)^2 \right)
+ \exp \left( -\frac{1}{2} \left( {z + h} \over \sigma_z \right)^2 \right) } \over { \sqrt{2\pi} \sigma_z } }
```

"""
function puff(scenario::Scenario, ::Type{GaussianPuff})

    stab = _stability(scenario)
    m = _release_mass(scenario)
    h = _release_height(scenario)
    u = _windspeed(scenario)

    return GaussianPuffSolution(
        scenario,  #scenario::Scenario
        :gaussian, #model::Symbol
        m,  #mass
        h,  #release height
        u,  #windspeed
        stab # stability class
    )

end


function (g::GaussianPuffSolution)(x,y,z,t)

    # domain check
    if (z<0)||(t<0)
        return 0.0
    end

    G = g.mass
    h = g.height
    u = g.windspeed
    stab = g.stability
    xc = abs(u*t) # location of center of cloud
    σx = downwind_dispersion(xc, Puff, stab)
    σy = crosswind_dispersion(xc, Puff, stab)
    σz = vertical_dispersion(xc, Puff, stab)

    c = ( G/((2*π)^(1.5)*σx*σy*σz)
        * exp(-0.5*((x-u*t)/σx)^2)
        * exp(-0.5*(y/σy)^2)
        * ( exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2) ) )

    c = isnan(c) ? 0.0 : c

    return c
end
