# defining type for dispatch
struct IntPuff <: PuffModel end

struct IntPuffSolution{F<:Number,N<:Number,S<:StabilityClass,E<:EquationSet} <: Puff
    scenario::Scenario
    model::Symbol
    rate::F
    mass_to_vol::F
    duration::F
    height::F
    windspeed::F
    npuffs::N
    stability::Type{S}
    equationset::E
end
IntPuffSolution(s,m,r,ρ,d,h,u,n,stab,es) = IntPuffSolution(s,m,promote(r,ρ,d,h,u,)...,n,stab,es)

@doc doc"""
    puff(::Scenario, IntPuff[, ::EquationSet]; kwargs...)

Returns the solution to an integrated Gaussian dispersion model, where the
release is modeled as a sequence of Gaussian puffs, for the given scenario.

```math
c\left(x,y,z,t\right) = \sum_{i}^{n-1} { {Q_{i,j} \Delta t} \over n }
{ { \exp \left( -\frac{1}{2} \left( {x - u \left( t - i \delta t \right) } \over \sigma_x \right)^2 \right) } \over { \sqrt{2\pi} \sigma_x } }
{ { \exp \left( -\frac{1}{2} \left( {y} \over \sigma_y \right)^2 \right) } \over { \sqrt{2\pi} \sigma_y } }\\
\times { { \exp \left( -\frac{1}{2} \left( {z - h} \over \sigma_z \right)^2 \right)
+ \exp \left( -\frac{1}{2} \left( {z + h} \over \sigma_z \right)^2 \right) } \over { \sqrt{2\pi} \sigma_z } }
```

where δt is Δt/n, and the σs are dispersion parameters correlated with the distance x. 
The `EquationSet` defines the set of correlations used to calculate the dispersion 
parameters and windspeed. The concentration returned is in volume fraction, assuming 
the puff is a gas at ambient conditions.

# Arguments
- `n::Integer`: the number of discrete gaussian puffs, defaults to infinity

"""
function puff(scenario::Scenario, ::Type{IntPuff}, eqs=DefaultPuffSet(); n::Number=Inf)

    stab = _stability(scenario)
    ṁ = _mass_rate(scenario)
    Δt = _duration(scenario)
    h = _release_height(scenario)
    u = windspeed(scenario,h,eqs)

    # jet at ambient conditions
    Tₐ = _atmosphere_temperature(scenario)
    Pₐ = _atmosphere_pressure(scenario)
    ρₐ = _gas_density(scenario.substance,Tₐ,Pₐ)

    if n == Inf
        return PalazziSolution(
            scenario,  #scenario::Scenario
            :intpuff, #model::Symbol
            plume(scenario, GaussianPlume, eqs),
            Δt,   # duration
            u,    # windspeed
            stab, # stability class
            eqs   # equation set
        )
    elseif n > 1
        return IntPuffSolution(
            scenario,  #scenario::Scenario
            :intpuff, #model::Symbol
            ṁ,    # massrate
            ρₐ,   # mass-to-vol
            Δt,   # duration
            h,    # release height
            u,    # windspeed
            n,    # number of puffs
            stab, # stability class
            eqs   # equation set
        )
    elseif n==1
        return GaussianPuffSolution(
            scenario,  #scenario::Scenario
            :gaussian, #model::Symbol
            ṁ*Δt, # mass
            ρₐ,   # mass_to_vol
            h,    # release height
            u,    # windspeed
            stab, # stability class
            eqs,  # equation set
        )
    else
        error("Number of puffs must be a positive integer value, or Inf")
    end
end


function (ip::IntPuffSolution{F,<:Integer,S,E})(x,y,z,t) where {F<:Number,S<:StabilityClass,E<:EquationSet}
    # domain check
    if (x<0)||(z<0)||(t<0)
        return zero(F)
    end

    Qi = ip.rate
    Δt = ip.duration
    n = ip.npuffs # number of intervals = number of puffs - 1
    h = ip.height
    u = ip.windspeed

    # Only account for puffs that have already been emitted
    Δt = min(t,Δt)

    # Volume per puff is the total volume divided by the number of puffs
    G = Qi*Δt/n

    # Gaussian dispersion in the x direction
    ∑g = 0
    δt = Δt/(n-1)
    for i in 0:(n-1)
        t′ = t-i*δt
        xc = u*t′ # center of cloud

        σx = downwind_dispersion(xc, S, ip.equationset)
        gx = t′>0 ? exp(-0.5*((x-u*t′)/σx)^2)/(√(2π)*σx) : 0

        σy = crosswind_dispersion(xc, S, ip.equationset)
        gy = exp(-0.5*(y/σy)^2)/(√(2π)*σy)

        σz = vertical_dispersion(xc, S, ip.equationset)
        gz = ( exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2) )/(√(2π)*σz)

        g = gx*gy*gz
        ∑g += isnan(g) ? zero(F) : g
    end

    # concentration
    c = G*∑g/ip.mass_to_vol

    return min(c,one(F))
end
