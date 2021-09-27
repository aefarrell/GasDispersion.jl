"""
    plume_rise(scenario::Scenario, plumerise)
Implements the Briggs plume rise equations for buoyancy and momentum driven
plume rise as described in the ISC3 model guide EPA-454/B-95-003b
"""
function plume_rise(scenario, plumerise)

    if plumerise==false
        return x->zero(x)
    end

    # physics parameters
    g = 9.80616 #m/s^2

    # parameters of the jet
    Dⱼ = scenario.jet_diameter
    uⱼ = scenario.jet_velocity
    Tᵣ = scenario.release_temperature

    # parameters of the environment
    u = scenario.windspeed
    Tₐ = scenario.ambient_temperature
    stability = scenario.pasquill_gifford

    if any(ismissing, [Dⱼ, uⱼ, Tᵣ, u, Tₐ, stability])
        error("Scenario is incomplete")
    end

    # buoyancy flux
    Fb = g * uⱼ * Dⱼ^2 * (Tᵣ - Tₐ) / (4Tᵣ)

    # momentum flux
    Fm = uⱼ^2 * Dⱼ^2 * Tₐ/(4Tᵣ)

    # stability check
    if stability ∈ Set(["A","B","C","D"])
        stable = false
    elseif stability == "E"
        stable = true
        Γ = 0.020    # default lapse rate K/m
        s = (g/Tₐ)*Γ # stability
    elseif stability == "F"
        stable = true
        Γ = 0.035    # default lapse rate K/m
        s = (g/Tₐ)*Γ # stability
    else
        err = string(stability_class, " is not a valid stability class")
        error(err)
    end

    if stable
        ΔTc = 0.019582*Tᵣ*uⱼ*√(s)
        if (Tᵣ - Tₐ) > ΔTc
            # buoyancy dominated plume rise
            final_rise = 2.6*(Fb/(u*s))^(1/3)
            xf = 2.0715*u/√(s)
            return x -> if (x < xf) min(1.60*(Fb*x^2/u^3)^(1/3), final_rise) else final_rise end
        else
            # momentum dominated plume rise
            stable_momentum_rise = 1.5*(Fm/(uⱼ*√(s)))^(1/3)
            unstable_momentum_rise = 3*Dⱼ*(uⱼ/u)
            final_rise = min(stable_momentum_rise, unstable_momentum_rise)
            xf = (π/2)*(u/√(s))
            β = (1/3) + (u/uⱼ)
            return x -> if (x < xf) min((3Fm*sin(x*√(s)/u)/(β^2*u*√(s)))^(1/3), final_rise) else final_rise end
        end
    else
        if Fb < 55
            ΔTc = 0.0297*Tᵣ*(uⱼ/Dⱼ)^(1/3)
            xf = 49*Fb^(5/8)
            buoyant_rise = 21.425*Fb^(3/4)/u
        else
            ΔTc = 0.00575*Tᵣ*(uⱼ^2/Dⱼ)^(1/3)
            xf = 119*Fb^(2/5)
            buoyant_rise = 38.71*Fb^(3/5)/u
        end

        if (Tᵣ - Tₐ) > ΔTc
            # buoyancy dominated plume rise
            final_rise = buoyant_rise
            return x -> if (x < xf) min(1.60*(Fb*x^2/u^3)^(1/3), final_rise) else final_rise end
        else
            # momentum dominated plume rise
            final_rise = 3Dⱼ*(uⱼ/u)
            xf = if (Fb<=0) 4Dⱼ*(uⱼ+3u)^2/(uⱼ*u) else xf end
            β = (1/3) + (u/uⱼ)
            return x -> if (x < xf) min((3Fm*x/(β*u)^2)^(1/3), final_rise) else final_rise end
        end
    end
end
