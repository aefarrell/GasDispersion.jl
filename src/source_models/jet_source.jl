struct JetSource <: SourceModel
    phase::Symbol
    cd::Number
    k::Number
    diameter::Number
    pressure::Number
    temperature::Number
    density::Number
    height::Number
    duration::Number
end
JetSource(; phase=:liquid, dischargecoef=0.63, k=1.4, diameter,
            pressure, temperature, density, height,
            duration=Inf) = JetSource(phase, dischargecoef, k, diameter,
            pressure, temperature, density, height, duration)


@doc doc"""
    scenario_builder(source::JetSource(kwargs), atmosphere::Atmosphere)
Returns returns a scenario with a simple jet source from a circular hole. The
jet can either be a liquid or a gas, which case it is assumed to be an ideal
gas and the jet is isentropic.

Liquid and gas discharge models are per *Guidelines for Consequence Analysis of
Chemical Release*, CCPS, New York (1999)

# Arguments
- `phase=:liquid`: the phase, either :liquid or :gas
- `dischargecoef::Number=0.61`: the discharge coefficient cd
- `k::Number=1.4`: the heat capacity ration cp/cv
- `diameter::Number`: the diameter of the hole
- `pressure::Number`: the pressure upstream of the jet
- `temperature::Number`: the temperature upstream of the jet
- `density::Number`: the density of the fluid upstream of the jet

"""
function scenario_builder( source::JetSource, atmosphere::Atmosphere)
    phase = source.phase
    duration = source.duration
    cd = source.cd
    k = source.k
    d  = source.diameter
    T₁ = source.temperature
    P₁ = source.pressure
    ρ₁  = source.density
    T₂ = atmosphere.temperature
    P₂ = atmosphere.pressure
    h  = source.height
    A  = (π/4)*d^2

    if phase==:liquid
        u = cd*√((2/ρ₁)*(P₁-P₂))
        m  = ρ₁*A*u
        uⱼ = u
        ρⱼ = ρ₁
        Pⱼ = P₂
        Tⱼ = T₁
    elseif phase==:gas
        # isentropic expansion, limited by choked flow
        η = max((P₂/P₁),(2/(k+1))^(k/(k-1)))
        ρu = cd*√(ρ₁*P₁*(2k/(k-1))*(η^(2/k) - η^((k+1)/k)))
        m = A*ρu
        ρⱼ = ρ₁*η^(1/k)
        uⱼ = ρu/ρⱼ
        Pⱼ = η*P₁
        Tⱼ = T₁*η^((k-1)/k)
    else
        err = "$phase is not a valid phase, try either :liquid or :gas"
        error(err)
    end

    r = Release(; mass_rate=m,
                  duration=duration,
                  diameter=d,
                  velocity=uⱼ,
                  height=h,
                  pressure=Pⱼ,
                  temperature=Tⱼ,
                  density=ρⱼ)
    return Scenario(r,atmosphere)
end
