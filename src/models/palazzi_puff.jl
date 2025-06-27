# defining type for dispatch
struct Palazzi <: PuffModel end

struct PalazziSolution{F<:Number,S<:StabilityClass,E<:EquationSet} <: Puff
    scenario::Scenario
    model::Symbol
    rate::F
    mass_to_vol::F
    duration::F
    height::F
    windspeed::F
    stability::Type{S}
    equationset::E
end
PalazziSolution(s,m,r,ρ,d,h,u,stab,es) = PalazziSolution(s,m,promote(r,ρ,d,h,u,)...,stab,es)

@doc doc"""
    puff(::Scenario, Palazzi[, ::EquationSet]; kwargs...)

...

# References
+ Palazzi, E, M De Faveri, Giuseppe Fumarola, and G Ferraiolo. “Diffusion from a Steady Source of Short Duration.” *Atmospheric Environment*. 16, no. 12 (1982): 2785–90.

"""
function puff(scenario::Scenario, ::Type{Palazzi}, eqs=DefaultPuffSet())

    stab = _stability(scenario)
    ṁ = _mass_rate(scenario)
    Δt = _duration(scenario)
    h = _release_height(scenario)
    u = windspeed(scenario,h,eqs)

    # jet at ambient conditions
    Tₐ = _atmosphere_temperature(scenario)
    Pₐ = _atmosphere_pressure(scenario)
    ρₐ = _gas_density(scenario.substance,Tₐ,Pₐ)

    return PalazziSolution(
        scenario,  #scenario::Scenario
        :intpuff, #model::Symbol
        ṁ,    # massrate
        ρₐ,   # mass-to-vol
        Δt,   # duration
        h,    # release height
        u,    # windspeed
        stab, # stability class
        eqs   # equation set
    )
end


function (ip::PalazziSolution{F,S,E})(x,y,z,t) where {F<:Number,S<:StabilityClass,E<:EquationSet}
    # domain check
    if (x<0)||(z<0)||(t<0)
        return zero(F)
    end

    Qi = ip.rate
    Δt = ip.duration
    h = ip.height
    u = ip.windspeed

    # Only account for puffs that have already been emitted
    Δt = min(t,Δt)

    # Gaussian dispersion in the x direction
    σx_a = downwind_dispersion(u*(t-Δt), S, ip.equationset)
    σx_b = downwind_dispersion(u*t, S, ip.equationset)
    a  = (x-u*(t-Δt))/(√2*σx_a)
    b  = (x-u*t)/(√2*σx_b)
    ∫gx = erf(b,a)/(2u)

    # Gaussian dispersion in the y direction
    σy = crosswind_dispersion(x, S, ip.equationset)
    gy = exp((-1/2)*(y/σy)^2)/(√(2π)*σy)

    # Gaussian dispersion in the z direction
    σz = vertical_dispersion(x, S, ip.equationset)
    gz = ( exp((-1/2)*((z-h)/σz)^2) + exp((-1/2)*((z+h)/σz)^2) )/(√(2π)*σz)

    # concentration
    c = Qi*∫gx*gy*gz/ip.mass_to_vol

    return min(c,one(F))
end