# for dispatch
struct GaussianMixingLayer <: PlumeModel end

# new vertical term
struct SimpleMixingLayer{N<:Integer,F<:Number} <: GaussianVerticalTerm
    nterms::N
    mixing_height::F
end

function vertical_term(z, h, σz, ml::SimpleMixingLayer) 
    Fz = exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2)
    for i=1:ml.nterms
        H₁ = z - (2*i*ml.mixing_height - h)
        H₂ = z + (2*i*ml.mixing_height - h)
        H₃ = z - (2*i*ml.mixing_height + h)
        H₄ = z + (2*i*ml.mixing_height + h)
        next_term = exp(-0.5*(H₁/σz)^2) + exp(-0.5*(H₂/σz)^2) + exp(-0.5*(H₃/σz)^2) + exp(-0.5*(H₄/σz)^2)
        Fz += next_term
        if next_term ≈ 0
            break
        end
    end
    return Fz/(√(2π)*σz)
end

@doc doc"""
    plume(::Scenario, GaussianMixingLayer[, ::EquationSet]; kwargs...)

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
+ 

"""
function plume(scenario::Scenario, ::Type{GaussianMixingLayer}, eqs=DefaultSet(); h_min=1.0, n_terms=10)
    # parameters of the jet
    hᵣ = _release_height(scenario)
    m  = _mass_rate(scenario)
    Q  = _release_flowrate(scenario)
    ρⱼ = _release_density(scenario)

    # jet at ambient conditions
    Tₐ = _atmosphere_temperature(scenario)
    Pₐ = _atmosphere_pressure(scenario)
    hₘ = _mixing_height(scenario)
    ρₐ = _gas_density(scenario.substance,Tₐ,Pₐ)
    Qᵒ = Q*ρⱼ/ρₐ

    # max concentration
    c_max = m/Qᵒ
    c_max = c_max/ρₐ

    return GaussianPlumeSolution(
    scenario,                               # scenario::Scenario
    :simplemixinglayer,                     # model::Symbol
    m,                                      # mass emission rate
    c_max,                                  # max_concentration
    ρₐ,                                     # mass concentration to vol concentration
    windspeed(scenario,max(hᵣ,h_min),eqs),  # windspeed
    hᵣ,                                     # effective_stack_height::Number
    SimpleCrossTerm(),
    SimpleMixingLayer(n_terms,hₘ),
    NoPlumeRise(),                          # plume rise model
    _stability(scenario),                   # stability class
    eqs                                     # equation set 
    )
end

@doc doc"""
    plume(::Scenario{AbstractSubstance,VerticalJet,Atmosphere}, GaussianPlume[, ::EquationSet]; kwargs...)

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
function plume(scenario::Scenario{<:AbstractSubstance,<:VerticalJet,<:Atmosphere}, ::Type{GaussianMixingLayer}, eqs=DefaultSet(); downwash::Bool=false, plumerise::Bool=true, h_min=1.0, n_terms=10)  
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
    hₘ = _mixing_height(scenario)

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
    :simplemixinglayer, #model::Symbol
    m,     #mass emission rate
    c_max, #max concentration
    ρₐ,    #mass to vol, at ambient conditions
    u,     #windspeed
    hᵣ,    #effective_stack_height::Number
    SimpleCrossTerm(),
    SimpleMixingLayer(n_terms,hₘ),
    plume, #plume rise model
    stab,  #stability class
    eqs    #equation set 
    )
end
