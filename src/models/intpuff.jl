# defining type for dispatch
struct IntPuff <: PuffModel end

struct IntPuffSolution{T<:Number,S<:StabilityClass} <: Puff
    scenario::Scenario
    model::Symbol
    massrate::Number
    duration::Number
    npuffs::T
    height::Number
    windspeed::Number
    stability::Type{S}
end

@doc doc"""
    puff(scenario::Scenario, IntPuff; kwargs...)

Generates an integrated gaussian dispersion model, where the release is modeled
as a sequence of gaussian puffs, for the given scenario and returns a
callable giving the concentration of the form `c(x, y, z, t)`

```math
c\left(x,y,z,t\right) = \sum_{i}^{n-1} { {\dot{m} \Delta t} \over n }
{ { \exp \left( -\frac{1}{2} \left( {x - u \left( t - i \delta t \right) } \over \sigma_x \right)^2 \right) } \over { \sqrt{2\pi} \sigma_x } }
{ { \exp \left( -\frac{1}{2} \left( {y} \over \sigma_y \right)^2 \right) } \over { \sqrt{2\pi} \sigma_y } }\\
\times { { \exp \left( -\frac{1}{2} \left( {z - h} \over \sigma_z \right)^2 \right)
+ \exp \left( -\frac{1}{2} \left( {z + h} \over \sigma_z \right)^2 \right) } \over { \sqrt{2\pi} \sigma_z } }
```

where δt is Δt/n, and the σs are dispersion parameters correlated with the distance x

# Arguments
- `n::Integer`: the number of discrete gaussian puffs, defaults to infinity

"""
function puff(scenario::Scenario, ::Type{IntPuff}; n::Number=Inf)

    stab = _stability(scenario)
    ṁ = _mass_rate(scenario)
    Δt = _duration(scenario)
    h = _release_height(scenario)
    u = _windspeed(scenario)

    if n > 1
        return IntPuffSolution(
            scenario,  #scenario::Scenario
            :intpuff, #model::Symbol
            ṁ,  #massrate
            Δt, #duration
            n, #number of puffs
            h,  #release height
            u,  #windspeed
            stab #stability class
        )
    elseif n==1
        return GaussianPuffSolution(
            scenario,  #scenario::Scenario
            :gaussian, #model::Symbol
            ṁ*Δt,  #mass
            h,  #release height
            u,  #windspeed
            stab #stability class
        )
    else
        error("Number of puffs must be a positive integer value, or Inf")
    end
end


function (ip::IntPuffSolution{<:Integer,<:StabilityClass})(x,y,z,t)
    # domain check
    if (z<0)||(t<0)
        return 0.0
    end

    ṁ = ip.massrate
    Δt = ip.duration
    n = ip.npuffs # number of intervals = number of puffs - 1
    h = ip.height
    u = ip.windspeed
    stab = ip.stability

    # Only account for puffs that have already been emitted
    Δt = min(t,Δt)

    # Mass per puff is the total mass divided by the number of puffs
    G = ṁ*Δt/n

    # Gaussian dispersion in the x direction
    ∑g = 0
    δt = Δt/(n-1)
    for i in 0:(n-1)
        t′ = t-i*δt
        xc = u*t′ # center of cloud

        σx = downwind_dispersion(xc, Puff, stab)
        gx = t′>0 ? exp(-0.5*((x-u*t′)/σx)^2)/(√(2π)*σx) : 0

        σy = crosswind_dispersion(xc, Puff, stab)
        gy = exp(-0.5*(y/σy)^2)/(√(2π)*σy)

        σz = vertical_dispersion(xc, Puff, stab)
        gz = ( exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2) )/(√(2π)*σz)

        g = gx*gy*gz
        ∑g += isnan(g) ? 0 : g
    end

    # Gaussian dispersion in the y and z directions


    # concentration
    c = G*∑g

    return c
end

function (ip::IntPuffSolution{Float64,<:StabilityClass})(x,y,z,t)
    # domain check
    if (z<0)||(t<0)
        return 0.0
    end

    ṁ = ip.massrate
    Δt = ip.duration
    h = ip.height
    u = ip.windspeed
    stab = ip.stability
    xc = u*t # center of cloud

    # Only account for puffs that have already been emitted
    Δt = min(t,Δt)

    # Gaussian dispersion in the x direction
    σx = downwind_dispersion(xc, Puff, stab)
    a  = (x-u*(t-Δt))/(√2*σx)
    b  = (x-u*t)/(√2*σx)
    ∫gx = erf(b,a)/(2u)

    # Gaussian dispersion in the y direction
    σy = crosswind_dispersion(xc, Puff, stab)
    gy = exp((-1/2)*(y/σy)^2)/(√(2π)*σy)

    # Gaussian dispersion in the z direction
    σz = vertical_dispersion(xc, Puff, stab)
    gz = ( exp((-1/2)*((z-h)/σz)^2) + exp((-1/2)*((z+h)/σz)^2) )/(√(2π)*σz)

    # concentration
    c = ṁ*∫gx*gy*gz

    return c
end
