struct BritterMcQuaidPlume <: PlumeModel end

struct BritterMcQuaidPlumeSolution <: Plume
    scenario::Scenario
    model::Symbol
    c₀::Number # initial concentration, kg/m³
    T′::Number # temperature correction
    D::Number  # critical length
    lb::Number # plume dimension parameter, m
    itp::Interpolations.GriddedInterpolation
    xnf::Number # near field distance
    xff::Number # far field distance
    A::Number  # far field constant
end

"""
    plume(::Scenario, BritterMcQuaidPlume)

Returns the solution to the Britter-McQuaid continuous ground level release
model for the given scenario.

# References
+ Britter, R.E. and J. McQuaid, *Workbook on the Dispersion of Dense Gases* HSE Contract Research Report No. 17/1988, 1988
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)

"""
function plume(scenario::Scenario, ::Type{BritterMcQuaidPlume})

    Q = _release_flowrate(scenario)
    ṁ = _mass_rate(scenario)
    ρⱼ = _release_density(scenario)
    Tᵣ = _release_temperature(scenario)

    u₁₀ = _windspeed(scenario, 10.0)
    ρₐ = _atmosphere_density(scenario)
    Tₐ = _atmosphere_temperature(scenario)

    # initial concentration
    c₀ = ṁ/Q

    # relative density
    const g = 9.80616  # gravity, m/s^2
    gₒ = g * ((ρⱼ - ρₐ)/ ρₐ)

    # critical distance
    D = √(Q/u₁₀)

    # temperature correction
    T′ = Tᵣ/Tₐ

    # correlation parameter
    α = 0.2*log10( gₒ^2 * Q * u₁₀^-5 )

    # setting up correlation for constant α
    # calculates the points for the linear interpolation
    concs, βs = _bm_pl_c(α)

    # near-field location
    xnf = 30
    xmin = 10^(minimum(βs))
    if xnf < xmin
        # linear interpolation between the near-field correlation
        # and the main correlation
        βnf = log10(xnf)
        cnf = 306 / (xnf^2 + 306)
        βs = [ βnf; βs]
        concs = [ cnf; concs ]
    else
        # for α > 0.605 the near-field overlaps with the correlation
        # stop the near-field correlation early
        # (this can lead to discontinuities)
        xnf = xmin
    end

    # linear interpolation
    itp = interpolate((βs,), concs, Gridded(Linear()))

    # far field correlation
    # starts at last interpolation point and decays like x′^-2
    xff = 10^(maximum(βs))
    A = minimum(concs)*xff^2

    # plume dimension parameters
    lb  = Q*gₒ/(u₁₀^3)

    return BritterMcQuaidPlumeSolution(
        scenario, #scenario::Scenario
        :brittermcquaid, #model::Symbol
        c₀,    # initial concentration, kg/m³
        T′,    # temperature_correction
        D,     # critical length, m
        lb,    # length parameter, m
        itp,   # interpolation
        xnf,   # near-field distance
        xff,   # far-field distance
        A      # far-field constant
    )
end

function (b::BritterMcQuaidPlumeSolution)(x, y, z)

    # plume dimension parameters
    Lu  = 0.5*b.D + 2*b.lb
    Lh0 = b.D + 8*b.lb

    # domain check
    if (x < -Lu) || (z < 0)
        return 0.0
    end

    # plume dimensions
    if x < 0
        # Britter-McQuaid model does not give a profile for x<0
        # assuming crosswind length is Lh0
        Lh = Lh0
    else
        Lh = Lh0 + 2.5*cbrt(b.lb*x^2)
    end

    # within the crosswind extent of the plume
    if (y < -Lh) || (y > Lh)
        return 0.0
    end

    x′ = x/b.D
    if x′ ≤ b.xnf
        # use near-field correlation
        c′ = x′ > 0 ? 306 / ( 306 + x′^2 ) : 1.0
    elseif b.xnf < x′ < b.xff
        # use linear interpolation
        β = log10(x′)
        c′ = b.itp(β)
    else
        # use far-field correlation
        # where A is a function of α
        c′ = b.A/(x′^2)
    end

    # non-isothermal correction
    c′ = c′ / (c′ + (1 - c′)*b.T′)
    c = b.c₀*c′

    # vertical extent, from continuity
    Lv = (b.D^2)/(2*Lh*c′)
    if z ≤ Lv
        return c
    else
        return 0.0
    end
end
