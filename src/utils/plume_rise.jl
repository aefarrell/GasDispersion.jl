abstract type PlumeRise end
abstract type BriggsModel <: PlumeRise end

struct NoPlumeRise <: PlumeRise end

struct BuoyantPlume <: BriggsModel
    Fb::Number
    xf::Number
    u::Number
    final_rise::Number
end

Base.isapprox(a::BuoyantPlume, b::BuoyantPlume) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

struct MomentumPlume{S<:StabilityClass} <: BriggsModel
    Fm::Number
    xf::Number
    β::Number
    s::Union{Number,Nothing}
    u::Number
    final_rise::Number
    stab::Type{S}
end

Base.isapprox(a::MomentumPlume, b::MomentumPlume) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

"""
    plume_rise(Dⱼ, uⱼ, Tᵣ, u, Tₐ, StabilityClass)
Implements the Briggs plume rise equations for buoyancy and momentum driven
plume rise.

# References
+ Briggs, G.A. *Plume Rise* U.S. Atomic Energy Commission, Oak Ridge (1969)
+ EPA, *User's Guide for the Industrial Source Complex (ISC3) Dispersion Models, vol 2*, U.S. Environmental Protection Agency EPA-454/B-95-003b (1995)

"""
function plume_rise(Dⱼ,uⱼ,Tᵣ,u,Tₐ, stab::Union{Type{ClassA},Type{ClassB},Type{ClassC},Type{ClassD}})
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
        return MomentumPlume(Fm,xf,β,nothing,u,final_rise,stab)
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
        xf = 2.0715*u/√(s)
        Δhf = 2.6*(Fb/(u*s))^(1/3)
        Δhf_calm = 4*(Fb^0.25)/(s^0.375)
        final_rise = min(Δhf,Δhf_calm)
        return BuoyantPlume(Fb,xf,u,final_rise)
    else
        # momentum dominated plume rise
        Fm = uⱼ^2 * Dⱼ^2 * Tₐ/(4Tᵣ)
        xf = (π/2)*(u/√(s))
        β = (1/3) + (u/uⱼ)
        Δhf = 1.5*(Fm/(uⱼ*√(s)))^(1/3)
        Δhf_unstable = 3*Dⱼ*(uⱼ/u)
        final_rise = min(Δhf,Δhf_unstable)
        return MomentumPlume(Fm,xf,β,s,u,final_rise,ClassE)
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
        xf = 2.0715*u/√(s)
        Δhf = 2.6*(Fb/(u*s))^(1/3)
        Δhf_calm = 4*(Fb^0.25)/(s^0.375)
        final_rise = min(Δhf,Δhf_calm)
        return BuoyantPlume(Fb,xf,u,final_rise)
    else
        # momentum dominated plume rise
        Fm = uⱼ^2 * Dⱼ^2 * Tₐ/(4Tᵣ)
        xf = (π/2)*(u/√(s))
        β = (1/3) + (u/uⱼ)
        Δhf = 1.5*(Fm/(uⱼ*√(s)))^(1/3)
        Δhf_unstable = 3*Dⱼ*(uⱼ/u)
        final_rise = min(Δhf,Δhf_unstable)
        return MomentumPlume(Fm,xf,β,s,u,final_rise,ClassF)
    end
end

function plume_rise(x::Number, m::BuoyantPlume)
    if x < m.xf
        return min(1.60*(m.Fb*x^2/m.u^3)^(1/3), m.final_rise)
    else
        return m.final_rise
    end
end

function plume_rise(x::Number, m::MomentumPlume{<:Union{ClassA,ClassB,ClassC,ClassD}})
    if x < m.xf
        return min((3m.Fm*x/(m.β*m.u)^2)^(1/3), m.final_rise)
    else
        return m.final_rise
    end
end

function plume_rise(x::Number, m::MomentumPlume{<:Union{ClassE,ClassF}})
    if x < m.xf
        return min((3m.Fm*sin(x*√(m.s)/m.u)/(m.β^2*m.u*√(m.s)))^(1/3), m.final_rise)
    else
        return m.final_rise
    end
end
