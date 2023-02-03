struct BritterMcQuaidPuff <: PuffModel end

struct BritterMcQuaidPuffSolution <: Puff
    scenario::Scenario
    model::Symbol
    c₀::Number    # initial concentration, kg/m³
    T′::Number    # temperature_correction
    V₀::Number    # initial volume, m³
    gₒ::Number    # reduced gravity, m/s²
    u₁₀::Number   # referebce windspeed, m/s
    itp::Interpolations.GriddedInterpolation
    xnf::Number   # near-field distance
    xff::Number   # far-field distance
    A::Number     # far-field constant
end

"""
    puff(scenario::Scenario, BritterMcQuaidPuff)

Returns the solution to the Britter-McQuaid instantaneous ground level release
model for the given scenario.

# References
+ Britter, R.E. and J. McQuaid, *Workbook on the Dispersion of Dense Gases* HSE Contract Research Report No. 17/1988, 1988
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)

"""
function puff(scenario::Scenario, ::Type{BritterMcQuaidPuff})

    Q = _release_flowrate(scenario)
    ṁ = _mass_rate(scenario)
    t = _duration(scenario)
    ρⱼ = _release_density(scenario)
    Tᵣ = _release_temperature(scenario)

    u₁₀ = _windspeed(scenario, 10.0)
    ρₐ = _atmosphere_density(scenario)
    Tₐ = _atmosphere_temperature(scenario)

    # initial cloud dimensions
    V₀ = Q*t

    # initial concentration
    c₀ = ṁ/Q

    # temperature correction
    T′ = Tᵣ/Tₐ

    # relative density
    g = 9.80616  # gravity, m/s^2
    gₒ = g * ((ρⱼ - ρₐ)/ ρₐ)

    # correlation parameter
    α = 0.5*log10( gₒ * cbrt(V₀) / u₁₀^2 )

    # setting up correlation for constant α
    # calculates the points for the linear interpolation
    concs, βs = _bm_pf_c(α)

    # near-field location
    xnf = 10^(minimum(βs))

    # linear interpolation
    itp = interpolate((βs,), concs, Gridded(Linear()))

    # far field correlation
    # starts at last interpolation point and decays like x′^-2
    xff = 10^(maximum(βs))
    A = minimum(concs)*xff^2

    return BritterMcQuaidPuffSolution(
        scenario, #scenario::Scenario
        :brittermcquaid, #model::Symbol
        c₀,    # initial concentration, kg/m³
        T′,    # temperature_correction
        V₀,    # initial volume, m³
        gₒ,    # reduced gravity, m/s²
        u₁₀,   # referebce windspeed, m/s
        itp,   # interpolation::Extrapolation
        xnf,   # near-field distance
        xff,   # far-field distance
        A      # far-field constant
    )

end

function (pf::BritterMcQuaidPuffSolution)(x,y,z,t)

    R₀ = cbrt(3*pf.V₀/4π)
    xc = 0.4*pf.u₁₀*t
    R² = R₀^2 + 1.2*√(pf.gₒ*pf.V₀)*t
    r² = (x-xc)^2 + y^2

    #domain check
    if r² > R²
        return 0.0
    end

    x′ = x/cbrt(pf.V₀)
    if x′ ≤ pf.xnf
        # use near-field correlation
        c′ = x′ > 0 ? 3.24 / (3.24 + x′^2) : 1.0

        # don't drop below the first correlation concentration
        c′ = max(c′, maximum(pf.itp.coefs))
    elseif pf.xnf < x′ < pf.xff
        # use linear interpolation
        β = log10(x′)
        c′ = pf.itp(β)
    else
        # use far-field correlation
        # where A is a function of α
        c′ = pf.A/(x′^2)
    end

    # non-isothermal correction
    c′ = c′ / (c′ + (1 - c′)*pf.T′)
    c = pf.c₀*c′

    # vertical extent, from continuity
    H = pf.V₀/(c′*π*R²)
    if z ≤ H
        return c
    else
        return 0.0
    end
end
