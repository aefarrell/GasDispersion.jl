const OMEGA = 7.2921e-5 # Earth's angular velocity in rad/s

_mixing_height(s::Scenario, eqs::BasicEquationSet) = _mixing_height(s.atmosphere, eqs)

function _mixing_height(a::SimpleAtmosphere{F,S}, es::BasicEquationSet; ϕ=40*π/180) where {F,S<:Union{UnstableClass,NeutralClass}}
    f = 2OMEGA*sin(ϕ) # Coriolis parameter
    u_str = friction_velocity(a, es)
    return 0.3*u_str/f
end
_mixing_height(a::SimpleAtmosphere{F,S}, es::BasicEquationSet) where {F,S<:StableClass} = Inf