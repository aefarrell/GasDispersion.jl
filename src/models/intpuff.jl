
# gaussian puff model
struct IntPuff{T<:Number} <: PuffModel
    npuffs::T
end

function IntPuff(;npuffs=Inf)
    if isfinite(npuffs)
        n = Int(npuffs)
        return IntPuff{Integer}(n)
    else
        return IntPuff{Float64}(Inf)
    end
end

struct IntPuffSolution{T<:Number} <: Puff
    scenario::Scenario
    model::Symbol
    massrate::Number
    duration::Number
    npuffs::T
    height::Number
    windspeed::Number
    downwind_dispersion::Dispersion
    crosswind_dispersion::Dispersion
    vertical_dispersion::Dispersion
end




@doc doc"""
    puff(scenario::Scenario, IntPuff(kwargs...))

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
- `npuffs::Integer`: the number of discrete gaussian puffs, defaults to infinity

"""
function puff(scenario::Scenario, model::IntPuff)

    stability = scenario.atmosphere.stability
    ṁ = scenario.release.mass_rate
    Δt = scenario.release.duration
    n = model.npuffs
    h = scenario.release.height
    u = scenario.atmosphere.windspeed

    # Pasquill-Gifford dispersion
    σx = crosswind_dispersion(stability; plume=false)
    σy = σx
    σz = vertical_dispersion(stability; plume=false)

    if n > 1
        return IntPuffSolution(
            scenario,  #scenario::Scenario
            :intpuff, #model::Symbol
            ṁ,  #massrate
            Δt, #duration
            n, #number of puffs
            h,  #release height
            u,  #windspeed
            σx, #downwind_dispersion::Dispersion
            σy, #crosswind_dispersion::Dispersion
            σz  #vertical_dispersion::Dispersion
        )
    elseif n==1
        return GaussianPuffSolution(
            scenario,  #scenario::Scenario
            :gaussian, #model::Symbol
            ṁ*Δt,  #mass
            h,  #release height
            u,  #windspeed
            σx, #downwind_dispersion::Dispersion
            σy, #crosswind_dispersion::Dispersion
            σz  #vertical_dispersion::Dispersion
        )
    else
        error("Number of puffs must be a positive integer value, or Inf")
    end
end


function (ip::IntPuffSolution{<:Integer})(x,y,z,t)
    ṁ = ip.massrate
    Δt = ip.duration
    n = ip.npuffs # number of intervals = number of puffs - 1
    h = ip.height
    u = ip.windspeed
    σx = ip.downwind_dispersion(x)
    σy = ip.crosswind_dispersion(x)
    σz = ip.vertical_dispersion(x)

    # Only account for puffs that have already been emitted
    Δt = min(t,Δt)

    # Mass per puff is the total mass divided by the number of puffs
    G = ṁ*Δt/n

    # Gaussian dispersion in the x direction
    gx = 0
    δt = Δt/(n-1)
    for i in 0:(n-1)
        t′ = t-i*δt
        gx′ = t′>0 ? exp(-0.5*((x-u*t′)/σx)^2)/(√(2π)*σx) : 0
        gx += isnan(gx′) ? 0 : gx′
    end

    # Gaussian dispersion in the y and z directions
    gy = exp(-0.5*(y/σy)^2)/(√(2π)*σy)
    gz = ( exp(-0.5*((z-h)/σz)^2) + exp(-0.5*((z+h)/σz)^2) )/(√(2π)*σz)

    # concentration
    c = G*gx*gy*gz

    return c
end

function (ip::IntPuffSolution{Float64})(x,y,z,t)
    ṁ = ip.massrate
    Δt = ip.duration
    h = ip.height
    u = ip.windspeed
    σx = ip.downwind_dispersion(x)
    σy = ip.crosswind_dispersion(x)
    σz = ip.vertical_dispersion(x)

    # Only account for puffs that have already been emitted
    Δt = min(t,Δt)

    # Gaussian dispersion in the x direction
    a  = (x-u*(t-Δt))/(√2*σx)
    b  = (x-u*t)/(√2*σx)
    ∫gx = erf(b,a)/(2u)

    # Gaussian dispersion in the y and z directions
    gy = exp((-1/2)*(y/σy)^2)/(√(2π)*σy)
    gz = ( exp((-1/2)*((z-h)/σz)^2) + exp((-1/2)*((z+h)/σz)^2) )/(√(2π)*σz)

    # concentration
    c = ṁ*∫gx*gy*gz

    return c
end
