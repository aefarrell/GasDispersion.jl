ambient = Dict{Symbol,Union{Missing, Number, String}}([
    :release_pressure => 101325.0,   # Pa
    :release_temperature => 298.15,# K
    :ambient_density => 1.225,     # kg/m^3
    :ambient_pressure => 101325.0,   # Pa
    :ambient_temperature => 298.15,# K
])

default_windspeeds = Dict{String, Number}([
# these are pretty arbitrary, it would be better to have a reference for these
    "A" => 3.0,
    "B" => 5.0,
    "C" => 6.0,
    "D" => 6.0,
    "E" => 5.0,
    "F" => 3.0,
])

"""
    scenario_builder(P, T, conditions; <keyword arguments>)

Builds a scenario given a release pressure, temperature, and set of ambient
conditions.

# Arguments
- P::Number, the release pressure in Pa
- T::Number, the release temperature in K
- conditions::Scenario=ambient, a scenario or dictionary with the ambient conditions
- model="jet", the model, defaults to a jet
- phase="liquid", the phase of the release, defaults to liquid
- stability::String, the Pasquill atmospheric stability
- windspeed::Number, the windspeed, defaults to a look-up table for stability
- kwargs... each model requires additional keyword arguments
"""
function scenario_builder(P::Number, T::Number, conditions::AbstractDict=ambient; stability::String, windspeed::Union{Missing,Number}=missing,
                          model::String="jet", phase::String="liquid", kwargs...)

    d = Dict{Symbol, Any}([ key => get(Missing, conditions, key) for key in fieldnames(Scenario)])

    if ismissing(windspeed)
        windspeed = default_windspeeds[stability]
    end

    d[:windspeed] = windspeed
    d[:pasquill_gifford] = stability

    if model == "jet" && phase == "liquid"
        required_params = [:liquid_density, :hole_diameter]
        if all(key -> key ∈ keys(kwargs), required_params)
            ρⱼ = kwargs[:liquid_density]
            Dⱼ = kwargs[:hole_diameter]
            cd = if (:discharge_coeff ∈ keys(kwargs)) kwargs[:discharge_coeff] else 0.61 end
            Pₐ = d[:ambient_pressure]

            m, uⱼ = liquid_jet(P, Pₐ, ρⱼ, Dⱼ; cd=cd)

            d[:mass_emission_rate] = m
            d[:jet_diameter] = Dⱼ
            d[:jet_velocity] = uⱼ
            d[:jet_density] = ρⱼ
            d[:release_pressure] = P
            d[:release_temperature] = T
            d[:release_height] = if (:release_height ∈ keys(kwargs)) kwargs[:release_height] else 0.0 end
        else
            missing_params = [ String(i) for i in filter(key -> !(key ∈ keys(kwargs)), required_params)]
            error_string = "Missing parameters: " * join(missing_params, ", ")
            error(error_string)
        end
    else
        error("Not implemented")
    end

    return Scenario(d)

end

function scenario_builder(P::Number, T::Number, conditions::Scenario; stability::String, windspeed::Union{Missing,Number}=missing,
                          model::String="jet", phase::String="liquid", kwargs...)

    d = Dict([ key => getproperty(ambient, key) for key in fieldnames(Scenario)])
    return scenario_builder(P, T; stability=stability, conditions=d, windspeed=windspeed,
                            model=model, phase=phase, kwargs...)
end


"""
    liquid_jet(P₁, P₂, ρ, D; cd=0.61)

Returns the mass emission rate (kg/s) and jet velocity (m/s) for a liquid jetting
through a circular hole of diameter D, with discharge coefficient cd.

# Arguments
- P₁::Number, the internal pressure in Pa
- P₂::Number, the external pressure in Pa
- ρ::Number, the liquid density in kg/m^3
- D::Number, the diameter of the hole in m, assumed circular
- cd::Number=0.61, the discharge coefficient
"""
function liquid_jet(P₁, P₂, ρ, D; cd=0.61)

    u = cd * √( 2*( (P₁ - P₂) / ρ))
    A = (π/4)*D^2
    m = ρ*A*u

    return m, u
end
