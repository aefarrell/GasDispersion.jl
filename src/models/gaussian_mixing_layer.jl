# for dispatch
struct GaussianMixingLayer <: PlumeModel end

# mixing layer vertical term -- method of images
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

# mixing layer vertical term -- periodic
struct PeriodicMixingLayer{N<:Integer,F<:Number} <: GaussianVerticalTerm
    nterms::N
    mixing_height::F
end

function vertical_term(z, h, σz, ml::PeriodicMixingLayer) 
    Fz = 1/2
    for n=1:ml.nterms
        next_term = cos(n*π*z/ml.mixing_height)*cos(n*π*h/ml.mixing_height)*exp(-0.5*(n*π*σz/ml.mixing_height)^2)
        Fz += next_term
        if next_term ≈ 0
            break
        end
    end
    return 2*Fz/ml.mixing_height
end

@doc doc"""
    plume(::Scenario, GaussianMixingLayer[, ::EquationSet]; **kwargs...)

Returns the solution to a Gaussian plume dispersion model with a simple reflective mixing layer.

```math
c(x, y, z) = \frac{m_i}{u} \frac{1}{ \sqrt{2\pi} \sigma_y} \exp\left(-\frac{1}{2}\left(\frac{y}{\sigma_y}\right)^2\right) F_z
```

where $ F\_z $ is the vertical dispersion term, a function of the mixing height $ h\_m $, and other symbols are as defined for a Gaussian
plume model.

There are two methods for the mixing layer:
- `:simplemixinglayer` uses a simple method of images, a series of reflections off of the mixing height
- `:periodicmixinglayer` uses an infinite series of cosine terms to calculate the vertical dispersion, which is more accurate for large mixing heights

# Keyword Arguments
- `h_min=1.0`:  Minimum height for windspeed calculations.
- `method=:simplemixinglayer`: The method used for the mixing layer.
- `n_terms=10`: Number terms in the series calculation of $ F_z $.
- `mixing_limit=10_000.0`: Limit for the mixing height, in meters. Mixing heights above this are treated as infinite.

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
+ US EPA. 1995. *User's Guide for the Industrial Source Complex (ISC3) Dispersion Models EPA-454/B-95-003b, vol 2*. Research Triangle Park, NC: Office of Air Quality Planning and Standards
+ Seinfeld, John H. and Spyros N. Pandis. 2006. *Atmospheric Chemistry and Physics* 2nd Ed. New York: John Wiley and Sons.

"""
function plume(scenario::Scenario, ::Type{GaussianMixingLayer}, eqs=DefaultSet(); 
               method=:simplemixinglayer, 
               h_min=1.0, 
               n_terms=10,
               mixing_limit=10_000.0)
    # parameters of the jet
    hᵣ = _release_height(scenario)
    m  = _mass_rate(scenario)
    Q  = _release_flowrate(scenario)
    ρⱼ = _release_density(scenario)

    # jet at ambient conditions
    Tₐ = _atmosphere_temperature(scenario)
    Pₐ = _atmosphere_pressure(scenario)
    hₘ = _mixing_height(scenario, eqs)
    ρₐ = _gas_density(scenario.substance,Tₐ,Pₐ)
    Qᵒ = Q*ρⱼ/ρₐ
    if hᵣ>hₘ
        error("Release height $hᵣ exceeds mixing height $hₘ. This is not supported by the Gaussian mixing layer model.")
    end

    # max concentration
    c_max = m/Qᵒ
    c_max = c_max/ρₐ

    # mixing layer
    if hₘ > mixing_limit
        # infinite mixing height
        @warn "Mixing height $hₘ exceeds limit $mixing_limit, using infinite mixing height."
        ml = SimpleVerticalTerm()
    elseif method == :simplemixinglayer
        ml = SimpleMixingLayer(n_terms, hₘ)
    elseif method == :periodicmixinglayer
        ml = PeriodicMixingLayer(n_terms, hₘ)
    else
        error("Unknown mixing layer method: $method")
    end

    # domain
    dom = ProblemDomain(0.0, Inf, -Inf, Inf, 0.0, hₘ)

    return GaussianPlumeSolution(
    scenario,                               # scenario::Scenario
    :simplemixinglayer,                     # model::Symbol
    m,                                      # mass emission rate
    c_max,                                  # max_concentration
    ρₐ,                                     # mass concentration to vol concentration
    windspeed(scenario,max(hᵣ,h_min),eqs),  # windspeed
    hᵣ,                                     # effective_stack_height::Number
    SimpleCrossTerm(),
    ml,
    NoPlumeRise(),                          # plume rise model
    _stability(scenario),                   # stability class
    eqs,                                    # equation set 
    dom)                                    # problem domain
end

@doc doc"""
    plume(::Scenario{Substance,VerticalJet,Atmosphere}, GaussianMixingLayer[, ::EquationSet]; **kwargs...)

Returns the solution to a Gaussian plume dispersion model with a simple reflective mixing layer.

```math
c(x, y, z) = \frac{m_i}{u} \frac{1}{ \sqrt{2\pi} \sigma_y} \exp\left(-\frac{1}{2}\left(\frac{y}{\sigma_y}\right)^2\right) F_z
```
where $ F\_z $ is the vertical dispersion term, a function of the mixing height $ h\_m $, and other symbols are as defined for a Gaussian
plume model.

There are two methods for the mixing layer:
- `:simplemixinglayer` uses a simple method of images, a series of reflections off of the mixing height
- `:periodicmixinglayer` uses an infinite series of cosine terms to calculate the vertical dispersion, which is more accurate for large mixing heights

# Keyword Arguments
- `downwash::Bool=false`: Include stack-downwash effects if true.
- `plumerise::Bool=true`: Include plume-rise effects using Briggs' model if true.
- `h_min=1.0`: Minimum height, in meters, for windspeed calculations.
- `method=:simplemixinglayer`: The method used for the mixing layer.
- `n_terms=10`: Number terms in the series calculation of $ F_z $.
- `mixing_limit=10_000.0`: Limit for the mixing height, in meters. Mixing heights above this are treated as infinite.

# References
- AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
- Briggs, Gary A. 1969. *Plume Rise* Oak Ridge: U.S. Atomic Energy Commission
- US EPA. 1995. *User's Guide for the Industrial Source Complex (ISC3) Dispersion Models EPA-454/B-95-003b, vol 2*. Research Triangle Park, NC: Office of Air Quality Planning and Standards
- Seinfeld, John H. and Spyros N. Pandis. 2006. *Atmospheric Chemistry and Physics* 2nd Ed. New York: John Wiley and Sons.

"""
function plume(scenario::Scenario{<:AbstractSubstance,<:VerticalJet,<:Atmosphere}, 
               ::Type{GaussianMixingLayer}, eqs=DefaultSet(); 
               downwash::Bool=false, plumerise::Bool=true,
               method=:simplemixinglayer,
               h_min=1.0, 
               n_terms=10,
               mixing_limit=10_000.0)  
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
    hₘ = _mixing_height(scenario, eqs)
    if hᵣ>hₘ
        error("Release height $hᵣ exceeds mixing height $hₘ. This is not supported by the Gaussian mixing layer model.")
    end

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
        if plume.final_rise > hₘ-hᵣ
            @warn "Plume rise $(plume.final_rise) exceeds mixing height $hₘ, plume rise is limited to the mixing height."
            plume = _new_plume_rise(plume, hₘ)
        end
    else
        plume = NoPlumeRise()
    end

    # mixing layer
    if hₘ > mixing_limit
        # infinite mixing height
        @warn "Mixing height $hₘ exceeds limit $mixing_limit, using infinite mixing height."
        ml = SimpleVerticalTerm()
    elseif method == :simplemixinglayer
        ml = SimpleMixingLayer(n_terms, hₘ)
    elseif method == :periodicmixinglayer
        ml = PeriodicMixingLayer(n_terms, hₘ)
    else
        error("Unknown mixing layer method: $method")
    end

    # domain
    dom = ProblemDomain(0.0, Inf, -Inf, Inf, 0.0, hₘ)

    return GaussianPlumeSolution(
    scenario, #scenario::Scenario
    :simplemixinglayer, #model::Symbol
    m,     #mass emission rate
    c_max, #max concentration
    ρₐ,    #mass to vol, at ambient conditions
    u,     #windspeed
    hᵣ,    #effective_stack_height::Number
    SimpleCrossTerm(),
    ml,
    plume, #plume rise model
    stab,  #stability class
    eqs,   #equation set 
    dom)   #problem domain
end
