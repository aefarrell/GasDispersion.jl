# Plume crosswind dispersion correlations
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

# Puff crosswind dispersion correlations
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassA}) = 0.18*x^0.92
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassB}) = 0.14*x^0.92
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassC}) = 0.10*x^0.92
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassD}) = 0.06*x^0.92
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassE}) = 0.04*x^0.92
crosswind_dispersion(x, ::Type{Puff}, ::Type{ClassF}) = 0.02*x^0.89

# Puff downwind dispersion correlations
function downwind_dispersion(x, ::Type{Puff}, stab::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD},Type{ClassE},Type{ClassF}})
    return crosswind_dispersion(x, Puff, stab)
end

# Plume vertical dispersion correlations
function vertical_dispersion(x, ::Type{Plume}, ::Type{ClassA})
    δ = 107.7
    β = -1.7172
    γ = 0.2770
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
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassA}) = 0.60*x^0.75
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassB}) = 0.53*x^0.73
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassC}) = 0.34*x^0.71
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassD}) = 0.15*x^0.70
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassE}) = 0.10*x^0.65
vertical_dispersion(x, ::Type{Puff}, ::Type{ClassF}) = 0.05*x^0.61
