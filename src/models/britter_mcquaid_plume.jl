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

# Britter-McQuaid plume correlations
# as digitzed in TSCREEN
function _bm_pl_c_1(α)
    # corresponds to (Cm/C0) = 0.10
    if α ≤ -0.55
        return 1.75
    elseif -0.55 < α ≤ -0.14
        return 0.24*α + 1.88
    elseif -0.14 < α ≤ 1.0
        return -0.5*α + 1.78
    else
        @warn "α= $α is out of range"
        return -0.5*α + 1.78
    end
end

function _bm_pl_c_05(α)
    # corresponds to (Cm/C0) = 0.05
    if α ≤ -0.68
        return 1.92
    elseif -0.68 < α ≤ -0.29
        return 0.36*α + 2.16
    elseif -0.29 < α ≤ -0.18
        return 2.06
    elseif -0.18 < α ≤ 1.0
        return -0.56*α + 1.96
    else
        @warn "α= $α is out of range"
        return -0.56*α + 1.96
    end
end

function _bm_pl_c_02(α)
    # corresponds to (Cm/C0) = 0.02
    if α ≤ -0.69
        return 2.08
    elseif -0.69 < α ≤ -0.31
        return 0.45*α + 2.39
    elseif -0.31 < α ≤ -0.16
        return 2.25
    elseif -0.16 < α ≤ 1.0
        return -0.54*α + 2.16
    else
        @warn "α= $α is out of range"
        return -0.54*α + 2.16
    end
end

function _bm_pl_c_01(α)
    # corresponds to (Cm/C0) = 0.01
    if α ≤ -0.70
        return 2.25
    elseif -0.70 < α ≤ -0.29
        return 0.49*α + 2.59
    elseif -0.29 < α ≤ -0.20
        return 2.45
    elseif -0.20 < α ≤ 1.0
        return -0.52*α + 2.35
    else
        @warn "α= $α is out of range"
        return -0.52*α + 2.35
    end
end

function _bm_pl_c_005(α)
    # corresponds to (Cm/C0) = 0.005
    if α ≤ -0.67
        return 2.40
    elseif -0.67 < α ≤ -0.28
        return 0.59*α + 2.80
    elseif -0.28 < α ≤ -0.15
        return 2.63
    elseif -0.15 < α ≤ 1.0
        return -0.48*α + 2.56
    else
        @warn "α= $α is out of range"
        return -0.48*α + 2.56
    end
end

function _bm_pl_c_002(α)
    # corresponds to (Cm/C0) = 0.002
    if α ≤ -0.69
        return 2.60
    elseif -0.69 < α ≤ -0.25
        return 0.39*α + 2.87
    elseif -0.25 < α ≤ -0.13
        return 2.77
    elseif -0.13 < α ≤ 1.0
        return -0.50*α + 2.71
    else
        @warn "α= $α is out of range"
        return -0.50*α + 2.71
    end
end

function _bm_pl_c(α)
    concs = [0.10, 0.05, 0.02, 0.01, 0.005, 0.002]
    βs = [
            _bm_pl_c_1(α),   # (C_m/C_0) = 0.1
            _bm_pl_c_05(α),  # (C_m/C_0) = 0.05
            _bm_pl_c_02(α),  # (C_m/C_0) = 0.02
            _bm_pl_c_01(α),  # (C_m/C_0) = 0.01
            _bm_pl_c_005(α), # (C_m/C_0) = 0.005
            _bm_pl_c_002(α)  # (C_m/C_0) = 0.002
    ]
    return concs, βs
end

"""
    plume(scenario::Scenario, BritterMcQuaidPlume)

Returns the solution to a Britter-McQuaid dispersion model for the given
scenario.

Currently only implements the max concentration at a downwind distance x, the
other coordinates are ignored.

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
    g = 9.80616  # gravity, m/s^2
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
        cnf = 306*xnf^-2 / (1+ 306*xnf^-2)
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
        itp,   # interpolation::Extrapolation
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
    if x < -Lu
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

    Lv = (b.D^2)/Lh
    # within the extent of the plume
    if (y < -Lh) || (y > Lh) || (z < 0) || (z > Lv)
        return 0.0
    end

    x′ = x/b.D
    if x′ ≤ b.xnf
        # use near-field correlation
        c′ = x′ > 0 ? 306*x′^-2 / (1+ 306*x′^-2) : 1.0
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
    c = b.c₀*c′ / (c′ + (1 - c′)*b.T′)
    return c
end
