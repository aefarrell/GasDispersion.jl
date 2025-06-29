const ω = 7.2921e-5 # Earth's angular velocity in rad/s

_mixing_height(s::Scenario, eqs::BasicEquationSet) = _mixing_height(s.atmosphere, eqs)

function _mixing_height(a::SimpleAtmosphere, es::BasicEquationSet; ϕ=40*π/180)
    stab = _stability(a)
    if stab ∈ (ClassA, ClassB, ClassC, ClassD)
        f = 2ω*sin(ϕ) # Coriolis parameter
        u_str = friction_velocity(a, es)
        return 0.3*u_str/f
    else
        return Inf # infinite mixing height for stable conditions
    end
end