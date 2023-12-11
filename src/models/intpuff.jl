# defining type for dispatch
struct IntPuff <: PuffModel end

struct IntPuffSolution{T<:Number,S<:StabilityClass} <: Puff
    scenario::Scenario
    model::Symbol
    rate::Number
    duration::Number
    npuffs::T
    height::Number
    windspeed::Number
    stability::Type{S}
    equationset::EquationSet
end

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

where δt is Δt/n, and the σs are dispersion parameters correlated with the distance x. The `EquationSet` defines the set of correlations used to calculate the dispersion parameters.

# Arguments
- `n::Integer`: the number of discrete gaussian puffs, defaults to infinity

"""
function puff(scenario::Scenario, ::Type{IntPuff}, eqs::EquationSet=DefaultSet(); n::Number=Inf)

    stab = _stability(scenario)
    ṁ = _mass_rate(scenario)
    ρⱼ = _release_density(scenario)
    Qi = ṁ/ρⱼ
    Δt = _duration(scenario)
    h = _release_height(scenario)
    u = _windspeed(scenario)

    if n > 1
        return IntPuffSolution(
            scenario,  #scenario::Scenario
            :intpuff, #model::Symbol
            Qi,    # massrate
            Δt,   # duration
            n,    # number of puffs
            h,    # release height
            u,    # windspeed
            stab, # stability class
            eqs   # equation set
        )
    elseif n==1
        return GaussianPuffSolution(
            scenario,  #scenario::Scenario
            :gaussian, #model::Symbol
            Qi*Δt, # mass
            h,    # release height
            u,    # windspeed
            stab, # stability class
            eqs,  # equation set
        )
    else
        error("Number of puffs must be a positive integer value, or Inf")
    end
end


function (ip::IntPuffSolution{<:Integer,S})(x,y,z,t) where {S<:StabilityClass}
    # domain check
    if (x<0)||(z<0)||(t<0)
        return 0.0
    end

    Qi = ip.rate
    Δt = ip.duration
    n = ip.npuffs # number of intervals = number of puffs - 1
    h = ip.height
    u = ip.windspeed
    eqs = ip.equationset

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

        σx = downwind_dispersion(xc, Puff, S, eqs)
        gx = t′>0 ? exp(-0.5*((x-u*t′)/σx)^2)/(√(2π)*σx) : 0

        σy = crosswind_dispersion(xc, Puff, S, eqs)
        gy = exp(-0.5*(y/σy)^2)/(√(2π)*σy)

        σz = vertical_dispersion(xc, Puff, S, eqs)
        gz = ( exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2) )/(√(2π)*σz)

        g = gx*gy*gz
        ∑g += isnan(g) ? 0 : g
    end

    # Gaussian dispersion in the y and z directions


    # concentration
    c = G*∑g

    return c
end

function (ip::IntPuffSolution{Float64,S})(x,y,z,t) where {S<:StabilityClass}
    # domain check
    if (x<0)||(z<0)||(t<0)
        return 0.0
    end

    Qi = ip.rate
    Δt = ip.duration
    h = ip.height
    u = ip.windspeed
    eqs = ip.equationset

    # Only account for puffs that have already been emitted
    Δt = min(t,Δt)

    # Gaussian dispersion in the x direction
    σx_a = downwind_dispersion(u*(t-Δt), Puff, S, eqs)
    σx_b = downwind_dispersion(u*t, Puff, S, eqs)
    a  = (x-u*(t-Δt))/(√2*σx_a)
    b  = (x-u*t)/(√2*σx_b)
    ∫gx = erf(b,a)/(2u)

    # Gaussian dispersion in the y direction
    σy = crosswind_dispersion(x, Puff, S, eqs)
    gy = exp((-1/2)*(y/σy)^2)/(√(2π)*σy)

    # Gaussian dispersion in the z direction
    σz = vertical_dispersion(x, Puff, S, eqs)
    gz = ( exp((-1/2)*((z-h)/σz)^2) + exp((-1/2)*((z+h)/σz)^2) )/(√(2π)*σz)

    # concentration
    c = Qi*∫gx*gy*gz

    return c
end
