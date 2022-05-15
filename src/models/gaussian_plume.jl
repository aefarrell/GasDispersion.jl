include("plume_rise.jl")

# gaussian plume model
struct GaussianPlume <: PlumeModel
    downwash::Bool
    plumerise::Bool
end
GaussianPlume(;downwash=false, plumerise=false) = GaussianPlume(downwash,plumerise)

# Solution to the gaussian plume
struct GaussianPlumeSolution <: Plume
    scenario::Scenario
    model::Symbol
    mass_rate::Number
    windspeed::Number
    effective_stack_height::Number
    plume_rise::Function
    crosswind_dispersion::Dispersion
    vertical_dispersion::Dispersion
end

@doc doc"""
    plume(::Scenario; GaussianPlume(kwargs...))

Generates a gaussian dispersion model for the given scenario and returns a
callable giving the concentration of the form `c(x, y, z[, t])`

Gaussian plume model is per *Guidelines for Consequence Analysis of Chemical
Release*, CCPS, New York (1999)

```math
c\left(x,y,z\right) = {G \over 2 \pi \sigma_{y} \sigma_{z} u}
\exp \left[ -\frac{1}{2} \left( y \over \sigma_{y} \right)^2 \right] \\
\times \left\{ \exp \left[ -\frac{1}{2} \left( { z -h } \over \sigma_{z} \right)^2 \right]
+ \exp \left[ -\frac{1}{2} \left( { z + h } \over \sigma_{z} \right)^2 \right] \right\}
```

where the σs are dispersion parameters correlated with the distance x

# Arguments
- `downwash::Bool=false`: when true, includes stack-downwash effects
- `plumerise::Bool=false`: when true, includes plume-rise effects using Briggs' model
"""
function plume(scenario::Scenario, model::GaussianPlume)
    # parameters of the jet
    G  = scenario.release.mass_rate
    Dⱼ = scenario.release.diameter
    uⱼ = scenario.release.velocity
    hᵣ = scenario.release.height

    # parameters of the environment
    u = scenario.atmosphere.windspeed
    stability = scenario.atmosphere.stability

    # stack-tip downwash check
    if (model.downwash==true) && (uⱼ < 1.5*u)
        Δh_dw = 2*Dⱼ*( (uⱼ/u) - 1.5 )
    else
        Δh_dw = 0.0
    end

    hᵣ = hᵣ + Δh_dw

    # plume rise
    Δh = plume_rise(scenario, model.plumerise)

    # Pasquill-Gifford dispersion
    σy = crosswind_dispersion(stability)
    σz = vertical_dispersion(stability)

    return GaussianPlumeSolution(
    scenario, #scenario::Scenario
    :gaussian, #model::Symbol
    G,  #mass emission rate
    u,  #windspeed
    hᵣ, #effective_stack_height::Number
    Δh, #plume_rise::Function
    σy, #crosswind_dispersion::Function
    σz  #vertical_dispersion::Function
    )

end

function (g::GaussianPlumeSolution)(x, y, z, t=0)
    G = g.mass_rate
    u = g.windspeed
    hᵣ = g.effective_stack_height
    Δh = g.plume_rise(x)
    σy = g.crosswind_dispersion(x)
    σz = g.vertical_dispersion(x)
    hₑ  = hᵣ + Δh
    σyₑ = √( (Δh/3.5)^2 + σy^2 )
    σzₑ = √( (Δh/3.5)^2 + σz^2 )

    c = ( G/(2*π*u*σyₑ*σzₑ)
        * exp(-0.5*(y/σyₑ)^2)
        * ( exp(-0.5*((z-hₑ)/σzₑ)^2) + exp(-0.5*((z+hₑ)/σzₑ)^2) ) )

    return c

end
