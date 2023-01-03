abstract type PlumeRise end
abstract type BriggsModel <: PlumeRise end

struct NoPlumeRise <: PlumeRise end

struct BuoyantPlume <: BriggsModel
    Fb::Number
    xf::Number
    u::Number
    final_rise::Number
end

struct MomentumPlume{S<:StabilityClass} <: BriggsModel
    Fm::Number
    xf::Number
    β::Number
    s::Number
    final_rise::Number
    stab::Type{S}
end

"""
    plume_rise(Dⱼ, uⱼ, Tᵣ, u, Tₐ, ::Type{Union{ClassA, ClassB, ClassC, ClassD}})
Implements the Briggs plume rise equations for buoyancy and momentum driven
plume rise as described in the ISC3 model guide EPA-454/B-95-003b
"""
function plume_rise(Dⱼ,uⱼ,Tᵣ,u,Tₐ, stab::Type{Union{ClassA, ClassB, ClassC, ClassD}})
    # physics parameters
    g = 9.80616 #m/s^2

    # buoyancy flux
    Fb = g * uⱼ * Dⱼ^2 * (Tᵣ - Tₐ) / (4Tᵣ)

    # momentum flux
    Fm = uⱼ^2 * Dⱼ^2 * Tₐ/(4Tᵣ)

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
        return BuoyantPlume(Fb,xf,u,buoyant_rise)
    else
        # momentum dominated plume rise
        # momentum flux
        Fm = uⱼ^2 * Dⱼ^2 * Tₐ/(4Tᵣ)
        xf = if (Fb<=0) 4Dⱼ*(uⱼ+3u)^2/(uⱼ*u) else xf end
        β = (1/3) + (u/uⱼ)
        final_rise = 3Dⱼ*(uⱼ/u)
        return MomentumPlume(Fm,xf,β,NaN,final_rise,stab)
    end
end

function plume_rise(Dⱼ,uⱼ,Tᵣ,u,Tₐ, ::Type{ClassE})
    # physics parameters
    g = 9.80616 #m/s^2

    Γ = 0.020 # default lapse rate K/m
    s = (g/Tₐ)*Γ # stability
    Fb = g * uⱼ * Dⱼ^2 * (Tᵣ - Tₐ) / (4Tᵣ) # buoyancy flux
    ΔTc = 0.019582*Tᵣ*uⱼ*√(s) # temperature cross

    if (Tᵣ - Tₐ) > ΔTc
        # buoyancy dominated plume rise
        final_rise = 2.6*(Fb/(u*s))^(1/3)
        xf = 2.0715*u/√(s)
        return BuoyantPlume(Fb,xf,u,final_rise)
    else
        # momentum dominated plume rise
        stable_momentum_rise = 1.5*(Fm/(uⱼ*√(s)))^(1/3)
        unstable_momentum_rise = 3*Dⱼ*(uⱼ/u)
        final_rise = min(stable_momentum_rise, unstable_momentum_rise)
        xf = (π/2)*(u/√(s))
        β = (1/3) + (u/uⱼ)
        return MomentumPlume(Fm,xf,β,s,final_rise,ClassE)
    end
end

function plume_rise(Dⱼ,uⱼ,Tᵣ,u,Tₐ, ::Type{ClassF})
    # physics parameters
    g = 9.80616 #m/s^2

    Γ = 0.035 # default lapse rate K/m
    s = (g/Tₐ)*Γ # stability
    Fb = g * uⱼ * Dⱼ^2 * (Tᵣ - Tₐ) / (4Tᵣ) # buoyancy flux
    ΔTc = 0.019582*Tᵣ*uⱼ*√(s) # temperature cross

    if (Tᵣ - Tₐ) > ΔTc
        # buoyancy dominated plume rise
        final_rise = 2.6*(Fb/(u*s))^(1/3)
        xf = 2.0715*u/√(s)
        return BuoyantPlume(Fb,xf,u,final_rise)
    else
        # momentum dominated plume rise
        stable_momentum_rise = 1.5*(Fm/(uⱼ*√(s)))^(1/3)
        unstable_momentum_rise = 3*Dⱼ*(uⱼ/u)
        final_rise = min(stable_momentum_rise, unstable_momentum_rise)
        xf = (π/2)*(u/√(s))
        β = (1/3) + (u/uⱼ)
        return MomentumPlume(Fm,xf,β,s,final_rise,ClassF)
    end
end

function plume_rise(x, m::BuoyantPlume)
    if x < m.xf
        return min(1.60*(m.Fb*x^2/m.u^3)^(1/3), m.final_rise)
    else
        return m.final_rise
    end
end

function plume_rise(x, m::MomentumPlume{<:Union{ClassA,ClassB,ClassC,ClassD}})
    if x < m.xf
        return min((3m.Fm*x/(m.β*u)^2)^(1/3), m.final_rise)
    else
        return m.final_rise
    end
end

function plume_rise(x, m::MomentumPlume{<:Union{ClassE,ClassF}})
    if x < m.xf
        return min((3m.Fm*sin(x*√(m.s)/m.u)/(m.β^2*m.u*√(m.s)))^(1/3), m.final_rise)
    else
        return m.final_rise
    end
end
