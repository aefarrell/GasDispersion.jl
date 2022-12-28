power_law_params = Dict(
    :A => (α=0.108),
    :B => (α=0.112),
    :C => (α=0.120),
    :D => (α=0.142),
    :E => (α=0.203),
    :F => (α=0.253)
)

"""
    _windspeed(scenario::Scenario, z::Number)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

# References

"""
function _windspeed(s::Scenario,z::Number)
    stab = _stability(s)
    if stab ∈ keys(power_law_params)
        α, = power_law_params[stab]
        u0 = _windspeed(s)
        z0 = _release_height(s)

        return u0*(z/z0)^α
    else
        err = "$stab is not a valid Pasquill-Gifford stability class"
        error(err)
    end
end


"""
    _windspeed(z, u, zR, λ, stability_class; k=0.35)
returns the windspeed function u(z) for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

# Arguments
`u` friction velocity
`zR` surface roughness
`λ`  Monin-Obukhov length
`stability_class` Pasquill stability class (A, B, C, D, E, F)
`k`  von Karman's constant, 0.35


# References
    Businger et al. 1971
"""
function _windspeed(z::Number, u::Number, zR::Number, λ::Number, stability_class::Symbol; k=0.35)
    if stability_class ∈ Set([:A,:B,:C])
        a(z) = (1-15*(z/λ))^0.25
        Ψ(a) = 2*log((1+a)/2) + log((1+a^2)/2) - 2*atan(a) + π/2
        return (u/k)*(log((z+zR)/zR) - Ψ(a(z)))
    elseif stability_class == :D
        return (u/k)*log((z+zR)/zR)
    elseif stability_class ∈ Set([:E, :F])
        return (u/k)*(log((z+zR)/zR) - 4.7*(z/λ))
    else
        err = "$stability_class is not a valid Pasquill-Gifford stability class"
        error(err)
    end
end
