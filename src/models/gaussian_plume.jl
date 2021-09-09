include("plume_rise.jl")

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
function gaussian_plume_factory(scenario; downwash=false, plumerise=false)
    # parameters of the jet
    Q = scenario.mass_emission_rate
    Dⱼ = scenario.jet_diameter
    uⱼ = scenario.jet_velocity
    hᵣ = scenario.release_height

    # parameters of the environment
    u = scenario.windspeed
    stability = scenario.pasquill_gifford

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


    function gaussian_plume(x, y, z, t=0)
        hₑ  = hᵣ + Δh(x)
        σyₑ = √( (Δh(x)/3.5)^2 + σy(x)^2 )
        σzₑ = √( (Δh(x)/3.5)^2 + σz(x)^2 )

        gaussian_plume = ( Q/(2*π*u*σyₑ*σzₑ)
                         * exp(-0.5*(y/σyₑ)^2)
                         *( exp(-0.5*((z-hₑ)/σzₑ)^2) + exp(-0.5*((z+hₑ)/σzₑ)^2) ))

    end

    return gaussian_plume

end
