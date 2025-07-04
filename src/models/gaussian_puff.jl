
# gaussian puff model
struct GaussianPuff <: PuffModel end

struct GaussianPuffSolution{F<:Number,E<:EquationSet} <: Puff
    scenario::Scenario
    model::Symbol
    mass::F
    mass_to_vol::F
    height::F
    windspeed::F
    equationset::E
end
GaussianPuffSolution(s,m,q,ρ,h,u,es) = GaussianPuffSolution(s,m,promote(q,ρ,h,u)...,es)

@doc doc"""
    puff(::Scenario, GaussianPuff[, ::EquationSet])

Returns the solution to a Gaussian puff dispersion model for the given scenario.

```math
c\left(x,y,z,t\right) = m_{i} \Delta t
{ { \exp \left( -\frac{1}{2} \left( {x - u t } \over \sigma_x \right)^2 \right) } \over { \sqrt{2\pi} \sigma_x } }
{ { \exp \left( -\frac{1}{2} \left( {y} \over \sigma_y \right)^2 \right) } \over { \sqrt{2\pi} \sigma_y } }\\
\times { { \exp \left( -\frac{1}{2} \left( {z - h} \over \sigma_z \right)^2 \right)
+ \exp \left( -\frac{1}{2} \left( {z + h} \over \sigma_z \right)^2 \right) } \over { \sqrt{2\pi} \sigma_z } }
```

where the σs are dispersion parameters correlated with the distance x. The 
`EquationSet` defines the set of correlations used to calculate the dispersion
parameters and windspeed. The concentration returned is in volume fraction, assuming the puff
is a gas at ambient conditions.

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
"""
function puff(scenario::Scenario, ::GaussianPuff, eqs=DefaultPuffSet(); h_min=1.0)
    m = _release_mass(scenario)
    Tₐ = _atmosphere_temperature(scenario)
    Pₐ = _atmosphere_pressure(scenario)
    ρₐ = _gas_density(scenario.substance,Tₐ,Pₐ)
    h = _release_height(scenario)
    u = windspeed(scenario,max(h,h_min),eqs)

    return GaussianPuffSolution(
        scenario,  #scenario::Scenario
        :gaussian, #model::Symbol
        m,    # mass
        ρₐ,   # mass-to-vol
        h,    # release height
        u,    # windspeed
        eqs   # equation set
    )

end


function (g::GaussianPuffSolution{F,E})(x,y,z,t) where {F,E}

    # domain check
    if (x<0)||(z<0)||(t<0)
        return zero(F)
    end

    G = g.mass
    h = g.height
    u = g.windspeed
    stab = _stability(g.scenario)
    xc = abs(u*t) # location of center of cloud
    σx = downwind_dispersion(xc, stab, g.equationset)
    σy = crosswind_dispersion(xc, stab, g.equationset)
    σz = vertical_dispersion(xc, stab, g.equationset)

    c = ( G/((2*π)^(1.5)*σx*σy*σz)
        * exp(-0.5*((x-u*t)/σx)^2)
        * exp(-0.5*(y/σy)^2)
        * ( exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2) ) )

    c_vol = isnan(c) ? zero(F) : c/g.mass_to_vol

    return min(c_vol,one(F))
end
