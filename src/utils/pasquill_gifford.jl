z_params = Dict(
    "A" => (δ=107.7, β=-1.7172, γ=0.2770),
    "B" => (δ=0.1355, β=0.8752, γ=0.0136),
    "C" => (δ=0.09623, β=0.9477, γ=-0.0020),
    "D" => (δ=0.04134, β=1.1737, γ=-0.0316),
    "E" => (δ=0.02275, β=1.3010, γ=-0.0450),
    "F" => (δ=0.01122, β=1.4024, γ=-0.0540)
)

y_params = Dict(
    "A" => (δ=0.423, β=0.9),
    "B" => (δ=0.313, β=0.9),
    "C" => (δ=0.210, β=0.9),
    "D" => (δ=0.136, β=0.9),
    "E" => (δ=0.102, β=0.9),
    "F" => (δ=0.0674, β=0.9)
)

"""
    crosswind_dispersion(stability_class::String)
returns the crosswind dispersion function σy(x) for a given Pasquill-Gifford
stability class
`x` is assumed to be in meters and `σy` is in meters
"""
function crosswind_dispersion(stability_class::String)
    if stability_class ∈ Set(["A","B","C","D","E","F"])
        δ, β = y_params[stability_class]
    else
        err = string(stability_class, " is not a valid stability class")
        error(err)
    end

    return x -> δ*x^β
end

"""
vertical_dispersion(stability_class::String)
returns the crosswind dispersion function σz(x) for a given Pasquill-Gifford
stability class
`x` is assumed to be in meters and `σz` is in meters
"""
function vertical_dispersion(stability_class::String)
    if stability_class ∈ Set(["A","B","C","D","E","F"])
        δ, β, γ = z_params[stability_class]
    else
        err = string(stability_class, " is not a valid stability class")
        error(err)
    end

    return x -> δ*(x^β)*exp(γ*log(x)^2)
end
