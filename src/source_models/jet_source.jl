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
    scenario_builder(source::JetSource(), atmosphere::Atmosphere=Ambient())
Returns a scenario
# Arguments
- P₁::Number, the internal pressure in Pa
- P₂::Number, the external pressure in Pa
- ρ::Number, the liquid density in kg/m^3
- D::Number, the diameter of the hole in m, assumed circular
- cd::Number=0.61, the discharge coefficient
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
