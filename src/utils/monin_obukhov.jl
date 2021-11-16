λ_params = Dict(
    "A" => (a=-11.4, b=0.10),
    "B" => (a=-26.0, b=0.17),
    "C" => (a=-123.0, b=0.30),
    #"D" => Monin-Obhukov length is infinite
    "E" => (a=123.0, b=0.30),
    "F" => (a=26.0, b=0.17)
)

"""
    monin_obukhov(stability::String, roughness::Number)
returns the Monin-Obukhov length for a given Pasquill-Gifford stability class
and surface roughness (in meters)
Curve fit from
    Pasquill, F., *Atmospheric Diffusion, 2nd Ed.*, Halstead Press, New York, 1974.
"""
function monin_obukhov(stability::String, roughness::Number)
    if stability ∈ Set(["A","B","C","E","F"])
        a, b = λ_params[stability]
        λ = a*roughness^b
    elseif stability == "D"
        λ = Inf
    else
        err = "$stability is not a valid Pasquill-Gifford stability class"
        error(err)
    end

    return λ
end
