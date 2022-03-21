include("plume_rise.jl")

struct GaussianPlume
    scenario::Scenario
    model::Symbol
    effective_stack_height::Number
    plume_rise::Function
    crosswind_dispersion::Dispersion
    vertical_dispersion::Dispersion
end

"""
    gaussian_plume_factory(scenario; downwash=false, plumerise=false)

Generates a gaussian dispersion model on the given scenario and returns a
function giving the concentration of the form
c(x, y, z[, t])

`downwash` controls whether or not stack-downwash is included, by default it
is not
`plumerise` controls whether or not plume height is adjusted, by default there
is no plume rise
"""
function gaussian_plume_factory(scenario::Scenario; downwash=false, plumerise=false)

    required_params = [:mass_emission_rate, :jet_diameter, :jet_velocity,
                       :release_height, :windspeed, :pasquill_gifford]
    if all(key -> !(ismissing(getproperty(scenario,key))), required_params)
        # parameters of the jet
        Q = scenario.mass_emission_rate
        Dⱼ = scenario.jet_diameter
        uⱼ = scenario.jet_velocity
        hᵣ = scenario.release_height

        # parameters of the environment
        u = scenario.windspeed
        stability = scenario.pasquill_gifford
    else
        missing_params = [ String(i) for i in filter(key -> ismissing(getproperty(scenario,key)), required_params)]
        error_string = "These parameters cannot be missing: " * join(missing_params, ", ")
        e = MissingException(error_string)
        throw(e)
    end

    # stack-tip downwash check
    if (downwash==true) && (uⱼ < 1.5*u)
        Δh_dw = 2*Dⱼ*( (uⱼ/u) - 1.5 )
    else
        Δh_dw = 0.0
    end

    hᵣ = hᵣ + Δh_dw

    # plume rise
    Δh = plume_rise(scenario, plumerise)

    # Pasquill-Gifford dispersion
    σy = crosswind_dispersion(stability)
    σz = vertical_dispersion(stability)

    return GaussianPlume(
    scenario, #scenario::Scenario
    :gaussian, #model::Symbol
    hᵣ, #effective_stack_height::Number
    Δh, #plume_rise::Function
    σy, #crosswind_dispersion::Function
    σz #vertical_dispersion::Function
    )

end

function (g::GaussianPlume)(x, y, z, t=0)
    Q = g.scenario.mass_emission_rate
    u = g.scenario.windspeed
    hᵣ = g.effective_stack_height
    Δh = g.plume_rise(x)
    σy = g.crosswind_dispersion(x)
    σz = g.vertical_dispersion(x)
    hₑ  = hᵣ + Δh
    σyₑ = √( (Δh/3.5)^2 + σy^2 )
    σzₑ = √( (Δh/3.5)^2 + σz^2 )

    gaussian_plume = ( Q/(2*π*u*σyₑ*σzₑ)
                     * exp(-0.5*(y/σyₑ)^2)
                     *( exp(-0.5*((z-hₑ)/σzₑ)^2) + exp(-0.5*((z+hₑ)/σzₑ)^2) ))

end
