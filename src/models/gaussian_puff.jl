
# gaussian puff model
struct GaussianPuff <: PuffModel end

struct GaussianPuffSolution{F<:Number,S<:StabilityClass,E<:EquationSet} <: Puff
    scenario::Scenario
    model::Symbol
    volume::F
    height::F
    windspeed::F
    stability::Type{S}
    equationset::Type{E}
end
GaussianPuffSolution(s,m,Q,h,u,stab,es) = GaussianPuffSolution(s,m,promote(Q,h,u)...,stab,es)

@doc doc"""
    puff(::Scenario, GaussianPuff[, ::EquationSet])

Returns the solution to a Gaussian puff dispersion model for the given scenario.

```math
c\left(x,y,z,t\right) = Q_{i,j} \Delta t
{ { \exp \left( -\frac{1}{2} \left( {x - u t } \over \sigma_x \right)^2 \right) } \over { \sqrt{2\pi} \sigma_x } }
{ { \exp \left( -\frac{1}{2} \left( {y} \over \sigma_y \right)^2 \right) } \over { \sqrt{2\pi} \sigma_y } }\\
\times { { \exp \left( -\frac{1}{2} \left( {z - h} \over \sigma_z \right)^2 \right)
+ \exp \left( -\frac{1}{2} \left( {z + h} \over \sigma_z \right)^2 \right) } \over { \sqrt{2\pi} \sigma_z } }
```

where the σs are dispersion parameters correlated with the distance x. The 
`EquationSet` defines the set of correlations used to calculate the dispersion
 parameters.

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
"""
function puff(scenario::Scenario, ::Type{GaussianPuff}, eqs=DefaultSet; h_min=1.0)

    stab = _stability(scenario)
    m = _release_mass(scenario)
    ρ = _release_density(scenario)
    V = m/ρ
    h = _release_height(scenario)
    u = _windspeed(scenario,max(h,h_min),eqs)

    return GaussianPuffSolution(
        scenario,  #scenario::Scenario
        :gaussian, #model::Symbol
        V,    # volume
        h,    # release height
        u,    # windspeed
        stab, # stability class
        eqs   # equation set
    )

end


function (g::GaussianPuffSolution{<:Number,S,E})(x,y,z,t) where {S<:StabilityClass,E<:EquationSet}

    # domain check
    if (x<0)||(z<0)||(t<0)
        return 0.0
    end

    G = g.volume
    h = g.height
    u = g.windspeed
    xc = abs(u*t) # location of center of cloud
    σx = downwind_dispersion(xc, Puff, S, E)
    σy = crosswind_dispersion(xc, Puff, S, E)
    σz = vertical_dispersion(xc, Puff, S, E)

    c = ( G/((2*π)^(1.5)*σx*σy*σz)
        * exp(-0.5*((x-u*t)/σx)^2)
        * exp(-0.5*(y/σy)^2)
        * ( exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2) ) )

    c = isnan(c) ? 0.0 : c

    return c
end
