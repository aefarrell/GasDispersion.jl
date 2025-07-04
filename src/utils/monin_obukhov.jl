function monin_obuknov(a::SimpleAtmosphere, es::BasicEquationSet)
    stab = _stability(a)
    wind = _wind_equation(es)
    zR = _surface_roughness(a)
    return monin_obuknov(zR,stab,wind)
end

"""
    monin_obuknov(roughness, StabilityClass, WindEquation)
returns the Monin-Obukhov length for a given Pasquill-Gifford stability class
and surface roughness (in meters)

# References
+ Pasquill, Frank. 1974. *Atmospheric Diffusion, 2nd Ed.* New York: Halstead Press, New York

"""
monin_obuknov(zR::Number, ::ClassA, ::Any) = -11.4*zR^0.10
monin_obuknov(zR::Number, ::ClassB, ::Any) = -26.0*zR^0.17
monin_obuknov(zR::Number, ::ClassC, ::Any) = -123.0*zR^0.30
monin_obuknov(zR::Number, ::ClassD, ::Any) = Inf
monin_obuknov(zR::Number, ::ClassE, ::Any) = 123.0*zR^0.30
monin_obuknov(zR::Number, ::ClassF, ::Any) = 26.0*zR^0.17

monin_obuknov(zR::Number, s::StabilityClass) = monin_obuknov(zR, s, nothing)