# for modularity and code re-use
abstract type GaussianCrossTerm end
abstract type GaussianVerticalTerm end

# for dispatch
struct GaussianPlume <: PlumeModel end

# Solution to the gaussian plume
struct GaussianPlumeSolution{F<:Number,C<:GaussianCrossTerm,V<:GaussianVerticalTerm,P<:PlumeRise,S<:StabilityClass,E<:EquationSet} <: Plume
    scenario::Scenario
    model::Symbol
    rate::F
    max_concentration::F
    mass_to_vol::F
    windspeed::F
    effective_stack_height::F
    crossterm::C
    verticalterm::V
    plumerise::P
    stability::Type{S}
    equationset::E
end
GaussianPlumeSolution(s,m,Q,c,ρ,u,h_eff,cross,vert,pr,stab,es) = GaussianPlumeSolution(s,m,promote(Q,c,ρ,u,h_eff)...,cross,vert,pr,stab,es)

struct SimpleCrossTerm <:  GaussianCrossTerm end
cross_term(y, σy, ::SimpleCrossTerm) = exp(-0.5*(y/σy)^2)/(√(2π)*σy)

struct SimpleVerticalTerm <: GaussianVerticalTerm end
vertical_term(z, h, σz, ::SimpleVerticalTerm) = ( exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2) )/(√(2π)*σz)


@doc doc"""
    plume(::Scenario, GaussianPlume[, ::EquationSet]; kwargs...)

Returns the solution to a Gaussian plume dispersion model for the given scenario.

```math
c\left(x,y,z\right) = {m_{i} \over { 2 \pi \sigma_{y} \sigma_{z} u } }
\exp \left[ -\frac{1}{2} \left( y \over \sigma_{y} \right)^2 \right] \\
\times \left\{ \exp \left[ -\frac{1}{2} \left( { z -h } \over \sigma_{z} \right)^2 \right]
+ \exp \left[ -\frac{1}{2} \left( { z + h } \over \sigma_{z} \right)^2 \right] \right\}
```

where the σs are dispersion parameters correlated with the distance x. The 
`EquationSet` defines the set of correlations used to calculate the dispersion 
parameters. The concentration returned is in volume fraction, assuming the plume 
is a gas at ambient conditions.

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers

"""
function plume(scenario::Scenario, ::Type{GaussianPlume}, eqs=DefaultSet(); h_min=1.0)
    # parameters of the jet
    hᵣ = _release_height(scenario)
    m  = _mass_rate(scenario)
    Q  = _release_flowrate(scenario)
    ρⱼ = _release_density(scenario)

    # jet at ambient conditions
    Tₐ = _atmosphere_temperature(scenario)
    Pₐ = _atmosphere_pressure(scenario)
    ρₐ = _gas_density(scenario.substance,Tₐ,Pₐ)
    Qᵒ = Q*ρⱼ/ρₐ

    # max concentration
    c_max = m/Qᵒ
    c_max = c_max/ρₐ

    return GaussianPlumeSolution(
    scenario,                               # scenario::Scenario
    :gaussian,                              # model::Symbol
    m,                                      # mass emission rate
    c_max,                                  # max_concentration
    ρₐ,                                      # mass concentration to vol concentration
    windspeed(scenario,max(hᵣ,h_min),eqs), # windspeed
    hᵣ,                                     # effective_stack_height::Number
    SimpleCrossTerm(),
    SimpleVerticalTerm(),
    NoPlumeRise(),                          # plume rise model
    _stability(scenario),                   # stability class
    eqs                                     # equation set 
    )
end

@doc doc"""
    plume(::Scenario{Substance,VerticalJet,Atmosphere}, GaussianPlume[, ::EquationSet]; kwargs...)

Returns the solution to a Gaussian plume dispersion model for a vertical jet. By default the Briggs
plume rise model is used.

```math
c\left(x,y,z\right) = {m_{i} \over { 2 \pi \sigma_{y} \sigma_{z} u } }
\exp \left[ -\frac{1}{2} \left( y \over \sigma_{y} \right)^2 \right] \\
\times \left\{ \exp \left[ -\frac{1}{2} \left( { z -h } \over \sigma_{z} \right)^2 \right]
+ \exp \left[ -\frac{1}{2} \left( { z + h } \over \sigma_{z} \right)^2 \right] \right\}
```

where the σs are dispersion parameters correlated with the distance x. The 
`EquationSet` defines the set of correlations used to calculate the dispersion 
parameters. The concentration returned is in volume fraction, assuming the plume 
is a gas at ambient conditions.

# Arguments
- `downwash::Bool=false`: when true, includes stack-downwash effects
- `plumerise::Bool=true`: when true, includes plume-rise effects using Briggs' model

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ Briggs, Gary A. 1969. *Plume Rise* Oak Ridge: U.S. Atomic Energy Commission

"""
function plume(scenario::Scenario{<:AbstractSubstance,<:VerticalJet,<:Atmosphere}, ::Type{GaussianPlume}, eqs=DefaultSet(); downwash::Bool=false, plumerise::Bool=true, h_min=1.0)  
    # parameters of the jet
    Dⱼ = _release_diameter(scenario)
    uⱼ = _release_velocity(scenario)
    hᵣ = _release_height(scenario)
    m  = _mass_rate(scenario)
    Q  = _release_flowrate(scenario)
    ρⱼ = _release_density(scenario)

    # parameters of the environment
    u = windspeed(scenario,max(hᵣ,h_min),eqs)
    stab = _stability(scenario)
    Γ = _lapse_rate(scenario)

    # jet at ambient conditions
    Tₐ = _atmosphere_temperature(scenario)
    Pₐ = _atmosphere_pressure(scenario)
    ρₐ = _gas_density(scenario.substance,Tₐ,Pₐ)
    Qᵒ = Q*ρⱼ/ρₐ

    # max concentration
    c_max = m/Qᵒ
    c_max = c_max/ρₐ

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
        plume = plume_rise(Dⱼ,uⱼ,Tᵣ,u,Tₐ,Γ,stab)
    else
        plume = NoPlumeRise()
    end

    return GaussianPlumeSolution(
    scenario, #scenario::Scenario
    :gaussian, #model::Symbol
    m,     #mass emission rate
    c_max, #max concentration
    ρₐ,    #mass to vol, at ambient conditions
    u,     #windspeed
    hᵣ,    #effective_stack_height::Number
    SimpleCrossTerm(),
    SimpleVerticalTerm(),
    plume, #plume rise model
    stab,  #stability class
    eqs    #equation set 
    )
end

function (g::GaussianPlumeSolution{F,C,V,NoPlumeRise,S,E})(x, y, z, t=0) where {F,C,V,S,E}
    # domain check
    h = g.effective_stack_height
    if (x==0)&&(y==0)&&(z==h)
        return g.max_concentration
    elseif (x≤0)||(z<0)
        return zero(F)
    else
        G = g.rate
        u = g.windspeed
        σy = crosswind_dispersion(x,S,g.equationset)
        σz = vertical_dispersion(x,S,g.equationset)

        Fy = cross_term(y,σy,g.crossterm)
        Fz = vertical_term(z,h,σz,g.verticalterm)
        c = (G/u)*Fy*Fz

        # c is in kg/m^3
        # use density at ambient conditions to convert to vol pct
        c_vol = c/g.mass_to_vol

        return min(c_vol,g.max_concentration)
    end
end

function (g::GaussianPlumeSolution{F,C,V,<:BriggsModel,S,E})(x, y, z, t=0) where {F,C,V,S,E}
    # domain check
    h = g.effective_stack_height
    if (x==0)&&(y==0)&&(z==h)
        return g.max_concentration
    elseif (x≤0)||(z<0)
        return zero(F)
    else
        G = g.rate
        u = g.windspeed
        h = g.effective_stack_height
        m = g.plumerise
        Δh = plume_rise(x, m)
        σy = crosswind_dispersion(x,S,g.equationset)
        σz = vertical_dispersion(x,S,g.equationset)
        hₑ  = h + Δh
        σyₑ = √( (Δh/3.5)^2 + σy^2 )
        σzₑ = √( (Δh/3.5)^2 + σz^2 )

        Fy = cross_term(y,σyₑ,g.crossterm)
        Fz = vertical_term(z,hₑ,σzₑ,g.verticalterm)
        c = (G/u)*Fy*Fz

        # c is in kg/m^3
        # use density at ambient conditions to convert to vol pct
        c_vol = c/g.mass_to_vol

        return min(c_vol,g.max_concentration)
    end
end
