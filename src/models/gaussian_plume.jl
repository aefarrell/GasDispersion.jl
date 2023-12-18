# for dispatch
struct GaussianPlume <: PlumeModel end

# Solution to the gaussian plume
struct GaussianPlumeSolution{F<:Number,P<:PlumeRise,S<:StabilityClass,E<:EquationSet} <: Plume
    scenario::Scenario
    model::Symbol
    max_concentration::F
    rate::F
    windspeed::F
    effective_stack_height::F
    plumerise::P
    stability::Type{S}
    equationset::Type{E}
end
GaussianPlumeSolution(s,m,c_max,Q,u,h_eff,pr,stab,es) = GaussianPlumeSolution(s,m,promote(c_max,Q,u,h_eff)...,pr,stab,es)

@doc doc"""
    plume(::Scenario, GaussianPlume[, ::EquationSet]; kwargs...)

Returns the solution to a Gaussian plume dispersion model for the given scenario.

```math
c\left(x,y,z\right) = { {Q_{i,j} \over { 2 \pi \sigma_{y} \sigma_{z} u } }
\exp \left[ -\frac{1}{2} \left( y \over \sigma_{y} \right)^2 \right] \\
\times \left\{ \exp \left[ -\frac{1}{2} \left( { z -h } \over \sigma_{z} \right)^2 \right]
+ \exp \left[ -\frac{1}{2} \left( { z + h } \over \sigma_{z} \right)^2 \right] \right\}
```

where the σs are dispersion parameters correlated with the distance x. The 
`EquationSet` defines the set of correlations used to calculate the dispersion 
parameters.

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers

"""
function plume(scenario::Scenario, ::Type{GaussianPlume}, eqs=DefaultSet; downwash::Bool=false, plumerise::Bool=false, h_min=1.0)
    # parameters of the jet
    ṁ  = _mass_rate(scenario)
    ρⱼ = _release_density(scenario)
    Dⱼ = _release_diameter(scenario)
    Qⱼ = _release_flowrate(scenario)
    uⱼ = _release_velocity(scenario)
    hᵣ = _release_height(scenario)

    # parameters of the environment
    u = _windspeed(scenario,max(hᵣ,h_min),eqs)
    stab = _stability(scenario)
    Γ = _lapse_rate(scenario)

    # max concentration
    Qi = ṁ/ρⱼ
    c_max = min(Qi/Qⱼ,1.0)

    return GaussianPlumeSolution(
    scenario, #scenario::Scenario
    :gaussian, #model::Symbol
    c_max, # max concentration
    Qi,     #mass emission rate
    u,     #windspeed
    hᵣ,    #effective_stack_height::Number
    NoPlumeRise(), #plume rise model
    stab,  #stability class
    eqs    #equation set 
    )
end

@doc doc"""
    plume(::Scenario{Substance,VerticalJet,Atmosphere}, GaussianPlume[, ::EquationSet]; kwargs...)

Returns the solution to a Gaussian plume dispersion model for the given scenario.

```math
c\left(x,y,z\right) = { {Q_{i,j} \over { 2 \pi \sigma_{y} \sigma_{z} u } }
\exp \left[ -\frac{1}{2} \left( y \over \sigma_{y} \right)^2 \right] \\
\times \left\{ \exp \left[ -\frac{1}{2} \left( { z -h } \over \sigma_{z} \right)^2 \right]
+ \exp \left[ -\frac{1}{2} \left( { z + h } \over \sigma_{z} \right)^2 \right] \right\}
```

where the σs are dispersion parameters correlated with the distance x. The 
`EquationSet` defines the set of correlations used to calculate the dispersion 
parameters.

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers

# Arguments
- `downwash::Bool=false`: when true, includes stack-downwash effects
- `plumerise::Bool=true`: when true, includes plume-rise effects using Briggs' model

"""
function plume(scenario::Scenario{<:Substance,<:VerticalJet,<:Atmosphere}, ::Type{GaussianPlume}, eqs=DefaultSet; downwash::Bool=false, plumerise::Bool=true, h_min=1.0)
    # parameters of the jet
    ṁ  = _mass_rate(scenario)
    ρⱼ = _release_density(scenario)
    Dⱼ = _release_diameter(scenario)
    Qⱼ = _release_flowrate(scenario)
    uⱼ = _release_velocity(scenario)
    hᵣ = _release_height(scenario)

    # parameters of the environment
    u = _windspeed(scenario,max(hᵣ,h_min),eqs)
    stab = _stability(scenario)
    Γ = _lapse_rate(scenario)

    # max concentration
    Qi = ṁ/ρⱼ
    c_max = min(Qi/Qⱼ,1.0)

    # stack-tip downwash check
    if (downwash==true) && (uⱼ < 1.5*u)
        Δh_dw = 2*Dⱼ*( (uⱼ/u) - 1.5 )
    else
        Δh_dw = 0.0
    end

    hᵣ = hᵣ + Δh_dw

    # plume rise
    if plumerise == true
        Tᵣ = _release_temperature(scenario)
        Tₐ = _atmosphere_temperature(scenario)
        plume = plume_rise(Dⱼ,uⱼ,Tᵣ,u,Tₐ,Γ,stab)
    else
        plume = NoPlumeRise()
    end

    return GaussianPlumeSolution(
    scenario, #scenario::Scenario
    :gaussian, #model::Symbol
    c_max, # max concentration
    Qi,     #mass emission rate
    u,     #windspeed
    hᵣ,    #effective_stack_height::Number
    plume, #plume rise model
    stab,  #stability class
    eqs    #equation set 
    )
end

function (g::GaussianPlumeSolution{<:Number,NoPlumeRise,S,E})(x, y, z, t=0) where {S<:StabilityClass,E<:EquationSet}
    # domain check
    h = g.effective_stack_height
    if (x==0)&&(y==0)&&(z==h)
        return g.max_concentration
    elseif (x≤0)||(z<0)
        return 0.0
    else
        G = g.rate
        u = g.windspeed
        σy = crosswind_dispersion(x,Plume,S,E)
        σz = vertical_dispersion(x,Plume,S,E)

        c = ( G/(2*π*u*σy*σz)
            * exp(-0.5*(y/σy)^2)
            * ( exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2) ) )

        return min(g.max_concentration,c)
    end
end

function (g::GaussianPlumeSolution{<:Number,<:BriggsModel,S,E})(x, y, z, t=0) where {S<:StabilityClass,E<:EquationSet}
    # domain check
    h = g.effective_stack_height
    if (x==0)&&(y==0)&&(z==h)
        return g.max_concentration
    elseif (x≤0)||(z<0)
        return 0.0
    else
        G = g.rate
        u = g.windspeed
        h = g.effective_stack_height
        m = g.plumerise
        Δh = plume_rise(x, m)
        σy = crosswind_dispersion(x,Plume,S,E)
        σz = vertical_dispersion(x,Plume,S,E)
        hₑ  = h + Δh
        σyₑ = √( (Δh/3.5)^2 + σy^2 )
        σzₑ = √( (Δh/3.5)^2 + σz^2 )

        c = ( G/(2*π*u*σyₑ*σzₑ)
            * exp(-0.5*(y/σyₑ)^2)
            * ( exp(-0.5*((z-hₑ)/σzₑ)^2) + exp(-0.5*((z+hₑ)/σzₑ)^2) ) )

        return min(g.max_concentration,c)
    end
end
