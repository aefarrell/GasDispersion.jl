struct JetSource <: SourceModel
    phase::Symbol
    cd::Number
    diameter::Number
    pressure::Number
    temperature::Number
    density::Number
    height::Number
end
JetSource(; phase=:liquid, dischargecoef=0.63, diameter,
            pressure, temperature, density,
            height) = JetSource(phase, dischargecoef, diameter,
            pressure, temperature, density, height)


"""
    scenario_builder(source::JetSource(kwargs), atmosphere::Atmosphere=Ambient())
Returns returns a scenario with a simple jet source from a circular hole. The
jet can either be a liquid or a gas, which case it is assumed to be an ideal
gas and the jet is isentropic.

# Arguments
- phase=:liquid         the phase of the release, either :liquid or :gas
- dischargecoef=0.63    the discharge coefficient Cd
- diameter              the diameter of the hole
- pressure              the pressure upstream of the jet
- temperature           the temperature upstream of the jet
- density               the density of the fluid upstream of the jet
"""
function scenario_builder( source::JetSource, atmosphere::Atmosphere=Ambient())
    cd = source.cd
    d  = source.diameter
    P₁ = source.pressure
    P₂ = atmosphere.pressure
    ρ  = source.density
    h  = source.height

    u = cd * √( 2*( (P₁ - P₂) / ρ))
    A = (π/4)*d^2
    m = ρ*A*u

    r = Release(; mass_rate=m,
                  duration=Inf,
                  diameter=d,
                  velocity=u,
                  height=h,
                  pressure=P₁,
                  temperature=source.temperature,
                  density=ρ)
    return Scenario(r,atmosphere)
end
