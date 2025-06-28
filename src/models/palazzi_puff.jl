# defining type for dispatch
struct Palazzi <: PuffModel end

struct PalazziSolution{F<:Number,S<:StabilityClass,E<:EquationSet} <: Puff
    scenario::Scenario
    model::Symbol
    disp::Symbol
    plume::Plume
    duration::F
    windspeed::F
    stability::Type{S}
    equationset::E
end

PalazziDefaultSet = BasicEquationSet{DefaultWind,Defaultσy,Defaultσy,Defaultσz}

downwind_dispersion(x::Number, stab::Any, ::Type{Defaultσy}) = crosswind_dispersion(x, stab, Defaultσy)

@doc doc"""
    puff(::Scenario, Palazzi[, ::EquationSet]; kwargs...)

Returns the solution to a short duration Gaussian puff dispersion model for the given scenario, based on the work of Palazzi *et al*.

```math
c\left(x,y,z,t\right) = \Chi\left(x,y,z\right) \frac{1}{2} \left( \mathrm{erf} \left( { {x - u (t-\Delta t)} \over \sqrt{2} \sigma_x } \right) - \mathrm{erf} \left( { {x - u t} \over \sqrt{2} \sigma_x } \right)  \right)
````

where Χ is a Gaussian plume model and the σs are dispersion parameters. The `EquationSet`
defines the set of correlations used to calculate the dispersion parameters and windspeed.
The concentration returned is in volume fraction, assuming the puff is a gas at ambient
conditions.

There are multiple variations of the Palazzi short duration model, changing how the
downwind dispersion, σx, is calculated:
- `:default` follows Palazzi and calculates σx at the downwind distance *x*
- `:intpuff` calculates σx at the downwind distance to the cloud center *u t*
- `:tno` follows the TNO Yellow Book eqn 4.60b, using the distance *x* while the plume is still attached and the distance to the cloud center, *ut*, afterwards

# Arguments
- `plume_model::Type{Plume} = GaussianPlume` : the plume model $\Chi$
- `disp_method = :default` : the method for calculating the downwind dispersion

# References
+ Palazzi, E, M De Faveri, Giuseppe Fumarola, and G Ferraiolo. “Diffusion from a Steady Source of Short Duration.” *Atmospheric Environment*. 16, no. 12 (1982): 2785–90.

"""
function puff(scenario::Scenario, ::Type{Palazzi}, eqs=PalazziDefaultSet(); plume_model=GaussianPlume, disp_method=:default)

    stab = _stability(scenario)
    Δt = _duration(scenario)
    h = _release_height(scenario)
    u = windspeed(scenario,h,eqs)

    return PalazziSolution(
        scenario,  #scenario::Scenario
        :palazzi, #model::Symbol
        disp_method, #dispersion model::Symbol
        plume(scenario, plume_model, eqs), # plume model
        Δt,   # duration
        u,    # windspeed
        stab, # stability class
        eqs   # equation set
    )
end


function (ps::PalazziSolution{F,S,E})(x,y,z,t) where {F<:Number,S<:StabilityClass,E<:EquationSet}
    # domain check
    if (x<0)||(z<0)||(t<0)
        return zero(F)
    end

    Χ = ps.plume(x,y,z,t)
    Δt = min(t,ps.duration)
    u = ps.windspeed

    if ps.disp == :default
        xa = xb = x
    elseif ps.disp == :intpuff
        xa = u*(t-Δt)
        xb = u*t
    elseif ps.disp == :tno
        if t < ps.duration
            xa = xb = x
        else
            xa = xb = u*t
        end
    end

    # Gaussian dispersion in the x direction
    σx_a = downwind_dispersion(xa, S, ps.equationset)
    σx_b = downwind_dispersion(xb, S, ps.equationset)
    a  = (x-u*(t-Δt))/(√2*σx_a)
     b  = (x-u*t)/(√2*σx_b)

    # concentration
    c = Χ*erf(b,a)/2

    return min(c,one(F))
end