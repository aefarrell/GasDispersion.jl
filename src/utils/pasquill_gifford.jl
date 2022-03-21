struct Dispersion
    direction::Symbol
    equation::Symbol
    δ::Number
    β::Number
    γ::Number
end

plume_z_params = Dict(
    "A" => (δ=107.7, β=-1.7172, γ=0.2770),
    "B" => (δ=0.1355, β=0.8752, γ=0.0136),
    "C" => (δ=0.09623, β=0.9477, γ=-0.0020),
    "D" => (δ=0.04134, β=1.1737, γ=-0.0316),
    "E" => (δ=0.02275, β=1.3010, γ=-0.0450),
    "F" => (δ=0.01122, β=1.4024, γ=-0.0540)
)

puff_z_params = Dict(
    "A" => (δ=0.60, β=0.75),
    "B" => (δ=0.53, β=0.73),
    "C" => (δ=0.34, β=0.71),
    "D" => (δ=0.15, β=0.70),
    "E" => (δ=0.10, β=0.65),
    "F" => (δ=0.05, β=0.61),
)

plume_y_params = Dict(
    "A" => (δ=0.423, β=0.9, tₐ=18.4),
    "B" => (δ=0.313, β=0.9, tₐ=18.4),
    "C" => (δ=0.210, β=0.9, tₐ=18.4),
    "D" => (δ=0.136, β=0.9, tₐ=18.3),
    "E" => (δ=0.102, β=0.9, tₐ=11.4),
    "F" => (δ=0.0674, β=0.9, tₐ=4.6)
)

puff_y_params = Dict(
    "A" => (δ=0.18, β=0.92),
    "B" => (δ=0.14, β=0.92),
    "C" => (δ=0.10, β=0.92),
    "D" => (δ=0.06, β=0.92),
    "E" => (δ=0.04, β=0.92),
    "F" => (δ=0.02, β=0.89),
)

"""
    crosswind_dispersion(stability_class::String; <keyword arguments>)
returns the crosswind dispersion function σy(x) for a given Pasquill-Gifford
stability class, `x` is assumed to be in meters and `σy` is in meters

# Arguments
`plume` determines if the dispersion is for a plume or instantaneous release
(puff), by default a plume is assumed.
`avg_time` the averaging time in seconds, falls back to an "instantaneous"
averaging time

# References
plume dispersion correlations are from
Spicer, T.O., and J.A. Havens, *User's Guide for the DEGADIS 2.1 Dense Gas
Dispersion Model*, EPA-450/4-89-019, November 1989, pp 45-46

puff dispersion correlations are from:
Slade, D.H., Diffusion from Instantaneous Sources. *Meteorology and Atomic
Energy*, TID-24190, USAEC, 1968, pp 163-175
"""
function crosswind_dispersion(stability_class::String; plume=true, avg_time=600)
    if stability_class ∈ Set(["A","B","C","D","E","F"])
        if plume
            δ, β, tₐ = plume_y_params[stability_class]
            δ = δ*(max(avg_time, tₐ)/600)^0.2
        else
            δ, β = puff_y_params[stability_class]
        end
    else
        err = "$stability_class is not a valid Pasquill-Gifford stability class"
        error(err)
    end

    return Dispersion(
        :crosswind, #direction::Symbol
        :eqn1,      #equation::Symbol
        δ,
        β,
        0.0
    )
end

"""
vertical_dispersion(stability_class::String, plume=true)
returns the vertical dispersion function σz(x) for a given Pasquill-Gifford
stability class
`x` is assumed to be in meters and `σz` is in meters
`plume` determines if the dispersion is for a plume or instantaneous release
(puff), by default a plume is assumed.

# References
plume dispersion correlations are from
Seinfeld, J.H., *Atmospheric Chemistry and Physics of Air Pollution*, John
Wiley and Sons, New York, 1986

puff dispersion correlations are from:
Slade, D.H., Diffusion from Instantaneous Sources. *Meteorology and
Atomic Energy*, TID-24190, USAEC, 1968, pp163-175
"""
function vertical_dispersion(stability_class::String; plume=true)
    if stability_class ∈ Set(["A","B","C","D","E","F"])
        if plume
            δ, β, γ = plume_z_params[stability_class]
            return Dispersion(
                :vertical, #direction::Symbol
                :eqn2,      #equation::Symbol
                δ,
                β,
                γ
            )
        else
            δ, β = puff_z_params[stability_class]
            return Dispersion(
                :vertical, #direction::Symbol
                :eqn1,      #equation::Symbol
                δ,
                β,
                0.0
            )
        end
    else
        err = "$stability_class is not a valid Pasquill-Gifford stability class"
        error(err)
    end
end


function(d::Dispersion)(x)
    if d.equation == :eqn1
        return d.δ*x^d.β
    elseif d.equation == :eqn2
        return d.δ*(x^d.β)*exp(d.γ*log(x)^2)
    else
        return 0
    end
end
