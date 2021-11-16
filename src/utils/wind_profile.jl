power_law_params = Dict(
    "A" => (α=0.108),
    "B" => (α=0.112),
    "C" => (α=0.120),
    "D" => (α=0.142),
    "E" => (α=0.203),
    "F" => (α=0.253)
)

"""
    windspeed(u0, z0, stability_class)
returns the windspeed function u(z) for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

# Arguments
`u0` reference windspeed
`z0` height of reference point
`stability_class` Pasquill stability class (A, B, C, D, E, F)

# References

"""
function windspeed(u0::Number, z0::Number, stability_class::String)
    if stability_class ∈ Set(["A","B","C","D","E","F"])
        α, = power_law_params[stability_class]
    else
        err = "$stability_class is not a valid Pasquill-Gifford stability class"
        error(err)
    end

    return z -> u0*(z/z0)^α
end


"""
    windspeed(u, zR, λ, stability_class; k=0.35)
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
function windspeed(u::Number, zR::Number, λ::Number, stability_class::String; k=0.35)
    if stability_class ∈ Set(["A","B","C"])
        a(z) = (1-15*(z/λ))^0.25
        Ψ(a) = 2*log((1+a)/2) + log((1+a^2)/2) - 2*atan(a) + π/2
        return z -> (u/k)*(log((z+zR)/zR) - Ψ(a(z)))
    elseif stability_class == "D"
        return z -> (u/k)*log((z+zR)/zR)
    elseif stability_class ∈ Set(["E", "F"])
        return z -> (u/k)*(log((z+zR)/zR) - 4.7*(z/λ))
    else
        err = "$stability_class is not a valid Pasquill-Gifford stability class"
        error(err)
    end
end
