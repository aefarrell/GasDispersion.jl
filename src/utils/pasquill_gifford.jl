# Plume crosswind dispersion functions
function crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassA}; avg_time=600.0)
    δ, β, tₐ = 0.423, 0.9, 18.4
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

function crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassB}; avg_time=600.0)
    δ, β, tₐ = 0.313, 0.9, 18.4
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

function crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassC}; avg_time=600.0)
    δ, β, tₐ = 0.210, 0.9, 18.4
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

function crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassD}; avg_time=600.0)
    δ, β, tₐ = 0.136, 0.9, 18.3
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

function crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassE}; avg_time=600.0)
    δ, β, tₐ = 0.102, 0.9, 11.4
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

function crosswind_dispersion(x, ::Type{Plume}, ::Type{ClassF}; avg_time=600.0)
    δ, β, tₐ = 0.0674, 0.9, 4.6
    δ = δ*(max(avg_time, tₐ)/600)^0.2
    return δ*x^β
end

# Puff crosswind dispersion functions
function crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassA})
    δ, β = 0.18, 0.92
    return δ*x^β
end

function crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassB})
    δ, β = 0.14, 0.92
    return δ*x^β
end

function crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassC})
    δ, β = 0.10, 0.92
    return δ*x^β
end

function crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassD})
    δ, β = 0.06, 0.92
    return δ*x^β
end

function crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassE})
    δ, β = 0.04, 0.92
    return δ*x^β
end

function crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassF})
    δ, β = 0.02, 0.89
    return δ*x^β
end

# Puff downwind dispersion functions
function downwind_dispersion(x, ::Type{Puff}, stab::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}})
    return crosswind_dispersion(x, Puff, stab)
end

# Plume vertical dispersion functions
function vertical_dispersion(x, ::Type{Plume}, ::Type{ClassA})
    δ, β, γ = 107.7, -1.7172, 0.2770
    return δ*(x^β)*exp(γ*log(x)^2)
end

function vertical_dispersion(x, ::Type{Plume}, ::Type{ClassB})
    δ=0.1355
    β=0.8752
    γ=0.0136
    return δ*(x^β)*exp(γ*log(x)^2)
end

function vertical_dispersion(x, ::Type{Plume}, ::Type{ClassC})
    δ=0.09623
    β=0.9477
    γ=-0.0020
    return δ*(x^β)*exp(γ*log(x)^2)
end

function vertical_dispersion(x, ::Type{Plume}, ::Type{ClassD})
    δ=0.04134
    β=1.1737
    γ=-0.0316
    return δ*(x^β)*exp(γ*log(x)^2)
end

function vertical_dispersion(x, ::Type{Plume}, ::Type{ClassE})
    δ=0.02275
    β=1.3010
    γ=-0.0450
    return δ*(x^β)*exp(γ*log(x)^2)
end

function vertical_dispersion(x, ::Type{Plume}, ::Type{ClassF})
    δ=0.01122
    β=1.4024
    γ=-0.0540
    return δ*(x^β)*exp(γ*log(x)^2)
end

# Puff vertical dispersion functions
function vertical_dispersion(x, ::Type{Puff}, ::Type{ClassA})
    δ=0.60
    β=0.75
    return δ*x^β
end

function vertical_dispersion(x, ::Type{Puff}, ::Type{ClassB})
    δ=0.53
    β=0.73
    return δ*x^β
end

function vertical_dispersion(x, ::Type{Puff}, ::Type{ClassC})
    δ=0.34
    β=0.71
    return δ*x^β
end

function vertical_dispersion(x, ::Type{Puff}, ::Type{ClassD})
    δ=0.15
    β=0.70
    return δ*x^β
end

function vertical_dispersion(x, ::Type{Puff}, ::Type{ClassE})
    δ=0.10
    β=0.65
    return δ*x^β
end

function vertical_dispersion(x, ::Type{Puff}, ::Type{ClassF})
    δ=0.05
    β=0.61
    return δ*x^β
end
