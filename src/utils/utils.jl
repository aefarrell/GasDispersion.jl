# methods to get properties
_temperature(a::Atmosphere) = a.T
_temperature(r::Release) = r.T
_atmosphere_temperature(s::Scenario) = _temperature(s.atmosphere)
_release_temperature(s::Scenario) = _temperature(s.release)

_pressure(a::Atmosphere) = a.P
_pressure(r::Release) = r.P
_atmosphere_pressure(s::Scenario) = _pressure(s.atmosphere)
_release_pressure(s::Scenario) = _pressure(s.release)

_mass_rate(r::Release) = r.ṁ
_mass_rate(s::Scenario) = _mass_rate(s.release)

_duration(r::Release) = r.Δt
_duration(s::Scenario) = _duration(s.release)

_mass(r::Release) = _mass_rate(r)*_duration(r)
_release_mass(s::Scenario) = _mass(s.release)

_diameter(r::Release) = r.d
_release_diameter(s::Scenario) = _diameter(s.release)

_area(r::Release) = (π/4)*r.d^2
_release_area(s::Scenario) = _area(s.release)

_velocity(a::Atmosphere) = a.u
_velocity(r::Release) = r.u
_release_velocity(s::Scenario) = _velocity(s.release)

_flowrate(r::Release) = _area(r)*_velocity(r)
_release_flowrate(s::Scenario) = _flowrate(s.release)

_height(a::Atmosphere) = a.h
_height(r::Release) = r.h
_release_height(s::Scenario) = _height(s.release)

_liquid_fraction(r::Release) = r.f_l
_release_liquid_fraction(s::Scenario) = _liquid_fraction(s.release)

_windspeed(a::Atmosphere) = a.u
_windspeed(s::Scenario) = _windspeed(s.atmosphere)
_windspeed(s::Scenario, z::Number, es::EquationSet=DefaultSet()) = _windspeed(s.atmosphere, z, es)
_windspeed_height(a::Atmosphere) = a.h
_windspeed_height(s::Scenario) = _windspeed_height(s.atmosphere)
_stability(a::Atmosphere) = a.stability
_stability(s::Scenario) = _stability(s.atmosphere)

# density functions
_liquid_density(s::Substance) = _liquid_density(s, s.T_ref, s.P_ref)
_liquid_density(s::Substance{<:Any,<:Number,<:Any,<:Any,<:Any}, T::Number, P::Number) = s.ρ_l
_liquid_density(s::Substance{<:Any,<:Function,<:Any,<:Any,<:Any}, T::Number, P::Number) = s.ρ_l(T,P)

_gas_density(s::Substance) = _gas_density(s, s.T_ref, s.P_ref)
_gas_density(s::Substance{<:Number,<:Any,<:Any,<:Any,<:Any}, T::Number, P::Number) = s.ρ_g*(s.T_ref/T)*(P/s.P_ref)
_gas_density(s::Substance{<:Function,<:Any,<:Any,<:Any,<:Any}, T::Number, P::Number) = s.ρ_g(T,P)

function _density(s::Substance, f_l, T, P)
    f_g = 1 - f_l
    ρ_l = _liquid_density(s,T,P)
    ρ_g = _gas_density(s,T,P)

    return 1/(f_l/ρ_l + f_g/ρ_g)
end

_density(a::DryAir, T, P) = P/(a.Rs*T)
_density(a::Atmosphere) = _density(a, _temperature(a), _pressure(a))
_atmosphere_density(s::Scenario) = _density(s.atmosphere)
_release_density(s::Scenario) = _density(s.substance, _release_liquid_fraction(s), _release_temperature(s), _release_pressure(s))

# correlations for Monin-Obukhov length
include("monin_obukhov.jl")

# correlations for windspeed
include("wind_profile.jl")

# Pasquill-Gifford dispersion correlations
include("pasquill_gifford.jl")

# Briggs' model for plume rise
include("plume_rise.jl")

# model correlations
include("britter_mcquaid_correls.jl")

## Equation Sets
# includes for defined sets of correlations go here
include("equation_sets/ccps.jl")
include("equation_sets/tno.jl")