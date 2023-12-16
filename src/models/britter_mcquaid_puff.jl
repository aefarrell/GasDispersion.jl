struct BritterMcQuaidPuff <: PuffModel end

struct BritterMcQuaidPuffSolution{F<:Number,I} <: Puff
    scenario::Scenario
    model::Symbol
    c₀::F    # initial concentration,
    T′::F    # temperature_correction
    V₀::F    # initial volume, m³
    gₒ::F    # reduced gravity, m/s²
    u₁₀::F   # reference windspeed, m/s
    xnf::F   # near-field distance
    xff::F   # far-field distance
    A::F     # far-field constant
    itp::I
end

"""
    puff(scenario::Scenario, BritterMcQuaidPuff[, equationset::EquationSet])

Returns the solution to the Britter-McQuaid instantaneous ground level release
model for the given scenario.

The `equationset` is used to calculate the windspeed at 10m, all other 
correlations are as per the Britter-McQuaid model. Unless otherwise specified
a default power-law wind profile is used.

# References
+ Britter, Rex E. and J. McQuaid. 1988. *Workbook on the Dispersion of Dense Gases. HSE Contract Research Report No. 17/1988*
+ AIChE/CCPS. 1999. *Guidelines for Consequence Analysis of Chemical Releases*. New York: American Institute of Chemical Engineers
"""
function puff(scenario::Scenario, ::Type{BritterMcQuaidPuff}, eqs=DefaultSet)

    Q = _release_flowrate(scenario)
    ṁ = _mass_rate(scenario)
    t = _duration(scenario)
    ρⱼ = _release_density(scenario)
    Tᵣ = _release_temperature(scenario)

    u₁₀ = _windspeed(scenario, 10.0, eqs)
    ρₐ = _atmosphere_density(scenario)
    Tₐ = _atmosphere_temperature(scenario)

    # initial cloud dimensions
    V₀ = Q*t

    # initial concentration
    Qi = ṁ/ρⱼ
    c₀ = min(Qi/Q,1.0)

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
    itp = LinearInterpolation(concs, βs)

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
        u₁₀,   # reference windspeed, m/s
        xnf,   # near-field distance
        xff,   # far-field distance
        A,     # far-field constant
        itp   # interpolation::Extrapolation
    )

end

function (pf::BritterMcQuaidPuffSolution)(x,y,z,t)

    R₀ = cbrt(3*pf.V₀/4π)
    xc = 0.4*pf.u₁₀*t
    R² = R₀^2 + 1.2*√(pf.gₒ*pf.V₀)*t
    r² = (x-xc)^2 + y^2

    #domain check
    if r² > R² || z < 0
        return 0.0
    end

    x′ = x/cbrt(pf.V₀)
    if x′ ≤ pf.xnf
        # use near-field correlation
        c′ = x′ > 0 ? 3.24 / (3.24 + x′^2) : 1.0

        # don't drop below the first correlation concentration
        c′ = max(c′, maximum(pf.itp.u))
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
