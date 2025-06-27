# defining type for dispatch
struct Palazzi <: PuffModel end

struct PalazziSolution{F<:Number,S<:StabilityClass,E<:EquationSet} <: Puff
    scenario::Scenario
    model::Symbol
    plume::Plume
    duration::F
    windspeed::F
    stability::Type{S}
    equationset::E
end
PalazziSolution(s,m,p,d,u,stab,es) = PalazziSolution(s,m,p,promote(d,u)...,stab,es)

@doc doc"""
    puff(::Scenario, Palazzi[, ::EquationSet]; kwargs...)

...

# References
+ Palazzi, E, M De Faveri, Giuseppe Fumarola, and G Ferraiolo. “Diffusion from a Steady Source of Short Duration.” *Atmospheric Environment*. 16, no. 12 (1982): 2785–90.

"""
function puff(scenario::Scenario, ::Type{Palazzi}, eqs=DefaultPuffSet(); plume_model=GaussianPlume)

    stab = _stability(scenario)
    Δt = _duration(scenario)
    h = _release_height(scenario)
    u = windspeed(scenario,h,eqs)

    return PalazziSolution(
        scenario,  #scenario::Scenario
        :palazzi, #model::Symbol
        plume(scenario, plume_model, eqs), # plume model
        Δt,   # duration
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

    Χ = ip.plume(x,y,z,t)

    Δt = ip.duration
    Δt = min(t,Δt)
    u = ip.windspeed

    # Gaussian dispersion in the x direction
    σx_a = downwind_dispersion(u*(t-Δt), S, ip.equationset)
    σx_b = downwind_dispersion(u*t, S, ip.equationset)
    a  = (x-u*(t-Δt))/(√2*σx_a)
    b  = (x-u*t)/(√2*σx_b)
    ∫gx = erf(b,a)/2

    # concentration
    c = Χ*∫gx

    return min(c,one(F))
end