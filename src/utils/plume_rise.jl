abstract type PlumeRise end
abstract type BriggsModel <: PlumeRise end

struct NoPlumeRise <: PlumeRise end

struct BuoyantPlume{F<:Number} <: BriggsModel
    Fb::F
    xf::F
    u::F
    final_rise::F
end
BuoyantPlume(Fb,xf,u,final) = BuoyantPlume(promote(Fb,xf,u,final)...)

Base.isapprox(a::BuoyantPlume, b::BuoyantPlume) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

struct MomentumPlume{F<:Number,S<:StabilityClass} <: BriggsModel
    Fm::F
    xf::F
    β::F
    u::F
    s::F
    final_rise::F
    stab::S
end
MomentumPlume(Fm,xf,β,u,final,s,stab) = MomentumPlume(promote(Fm,xf,β,u,final,s)...,stab)

Base.isapprox(a::MomentumPlume, b::MomentumPlume) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])


# some helper functions to limit the plume rise, e.g. when a mixing layer prevents the plume from rising too high
"""
    _new_plume_rise(plume::BriggsModel, max_rise::Number)
Creates a new plume rise model with the same parameters as `plume`, but limits the final
rise to `max_rise`. This is useful when the plume rise exceeds the mixing height or
other constraints in the model.
"""
_new_plume_rise(plume::BuoyantPlume, max_rise::Number) = BuoyantPlume(plume.Fb, plume.xf, plume.u, max_rise)
_new_plume_rise(plume::MomentumPlume, max_rise::Number) = MomentumPlume(plume.Fm, plume.xf, plume.β, plume.u, plume.s, max_rise, plume.stab)

"""
    plume_rise(Dⱼ,uⱼ,Tᵣ,a::Atmosphere)
Implements the Briggs plume rise equations for buoyancy and momentum driven
plume rise.

# References
+ Briggs, Gary A. 1969. *Plume Rise* Oak Ridge: U.S. Atomic Energy Commission
+ EPA. 1995. *User's Guide for the Industrial Source Complex (ISC3) Dispersion Models, vol 2*. United States Environmental Protection Agency EPA-454/B-95-003b

"""
function plume_rise(Dⱼ,uⱼ,Tᵣ,u,Tₐ,Γ,stab::Union{UnstableClass,NeutralClass})
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
        return MomentumPlume(Fm,xf,β,u,0.0,final_rise,stab)
    end
end

function plume_rise(Dⱼ,uⱼ,Tᵣ,u,Tₐ,Γ,stab::StableClass)
    # physics parameters
    g = 9.80616 #m/s^2
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
        return MomentumPlume(Fm,xf,β,u,s,final_rise,stab)
    end
end

function plume_rise(x::Number, m::BuoyantPlume)
    if x < m.xf
        return min(1.60*(m.Fb*x^2/m.u^3)^(1/3), m.final_rise)
    else
        return m.final_rise
    end
end

function plume_rise(x::Number, m::MomentumPlume{<:Number,<:Union{UnstableClass,NeutralClass}})
    if x < m.xf
        return min((3m.Fm*x/(m.β*m.u)^2)^(1/3), m.final_rise)
    else
        return m.final_rise
    end
end

function plume_rise(x::Number, m::MomentumPlume{<:Number,<:StableClass})
    if x < m.xf
        return min((3m.Fm*sin(x*√(m.s)/m.u)/(m.β^2*m.u*√(m.s)))^(1/3), m.final_rise)
    else
        return m.final_rise
    end
end
