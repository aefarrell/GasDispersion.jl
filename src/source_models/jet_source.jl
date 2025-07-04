# defining type for dispatch
struct JetSource <: SourceModel end

@doc doc"""
    scenario_builder(substance::AbstractSubstance, ::JetSource, atmosphere::Atmosphere; kwargs...)
Returns returns a scenario for a simple jet from a circular hole. The
jet can either be a liquid or a gas (in which case it is assumed to be an ideal
gas and the jet is isentropic).

# References
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers

# Arguments
- `phase=:liquid`: the phase, either :liquid or :gas
- `dischargecoef::Number=0.61`: the discharge coefficient cd
- `diameter::Number`: the diameter of the hole, m
- `height::Number`: the height of the hole, m
- `pressure::Number`: the pressure upstream of the jet, Pa
- `temperature::Number`: the temperature upstream of the jet, K
- `duration::Number`: the duration of the leak, s
- `jet=:horizontal`: the type of jet, either :horizontal or :vertical

"""
function scenario_builder(substance::AbstractSubstance, ::JetSource, atmosphere::Atmosphere;
                          phase=:liquid,dischargecoef=0.63,diameter,pressure,temperature,height,duration=Inf,
                          jet=:horizontal)
    cd = dischargecoef
    d  = diameter
    T₁ = temperature
    P₁ = pressure
    T₂ = _temperature(atmosphere)
    P₂ = _pressure(atmosphere)
    k = substance.k
    h  = height
    A  = (π/4)*d^2

    if (h < 0) error("Height must be a positive number, instead got h = $h") end

    if phase==:liquid
        ρ₁ = _liquid_density(substance, T₁, P₁)
        u = cd*√((2/ρ₁)*(P₁-P₂))
        m  = ρ₁*A*u
        uⱼ = u
        ρⱼ = ρ₁
        Pⱼ = P₂
        Tⱼ = T₁
        f_l = 1.0
    elseif phase==:gas
        # isentropic expansion, limited by choked flow
        ρ₁ = _gas_density(substance, T₁, P₁)
        η = max((P₂/P₁),(2/(k+1))^(k/(k-1)))
        ρu = cd*√(ρ₁*P₁*(2k/(k-1))*(η^(2/k) - η^((k+1)/k)))
        m = A*ρu
        ρⱼ = ρ₁*η^(1/k)
        uⱼ = ρu/ρⱼ
        Pⱼ = η*P₁
        Tⱼ = T₁*η^((k-1)/k)
        f_l = 0.0
    else
        error("$phase is not a valid phase, try either :liquid or :gas")
    end

    if jet == :horizontal
        r = HorizontalJet(; mass_rate=m,
                            duration=duration,
                            diameter=d,
                            velocity=uⱼ,
                            height=h,
                            pressure=Pⱼ,
                            temperature=Tⱼ,
                            fraction_liquid=f_l)
    elseif jet == :vertical
        r = VerticalJet(; mass_rate=m,
                          duration=duration,
                          diameter=d,
                          velocity=uⱼ,
                          height=h,
                          pressure=Pⱼ,
                          temperature=Tⱼ,
                          fraction_liquid=f_l)
    else
        error("$jet is not a valid release type, try either :horizontal or :vertical")
    end

    return Scenario(substance,r,atmosphere)
end
