struct BritterMcQuaidPlume <: PlumeModel end

struct BritterMcQuaidPlumeSolution <: Plume
    scenario::Scenario
    model::Symbol
    ρⱼ::Number # jet density
    T′::Number # temperature correction
    D::Number  # critical length
    lb::Number # plume dimension parameter, m
    itp::Interpolations.Extrapolation
end

Britter_McQuaid_correlations = Dict(
    0.010 => ( αs=[-1.0, -0.7, -0.29, -0.2, 1.0],
               βs=[2.25, 2.25, 2.45, 2.45, 1.83]),
    0.005 => ( αs=[-1.0, -0.67, -0.28, -0.15, 1.0],
               βs=[2.4, 2.4, 2.63, 2.63, 2.07]),
    0.020 => ( αs=[-1.0, -0.69, -0.31, -0.16, 1.0],
               βs=[2.08, 2.08, 2.25, 2.25, 1.62]),
    0.002 => ( αs=[-1.0, -0.69, -0.25, -0.13, 1.0],
               βs=[2.6, 2.6, 2.77, 2.77, 2.21]),
    0.100 => ( αs=[-1.0, -0.55, -0.14, 1.0],
               βs=[1.75, 1.75, 1.85, 1.28]),
    0.050 => ( αs=[-1.0, -0.68, -0.29, -0.18, 1.0],
               βs=[1.92, 1.92, 2.06, 2.06, 1.4]),
)


"""
    plume(scenario::Scenario, BritterMcQuaidPlume)

Returns the solution to a Britter-McQuaid dispersion model for the given
scenario.

Currently only implements the max concentration at a downwind distance x, the
other coordinates are ignored.

"""
function plume(scenario::Scenario, ::Type{BritterMcQuaidPlume})

    Q = _mass_rate(scenario)
    h = _release_height(scenario)
    ρⱼ = _release_density(scenario)
    Tᵣ = _release_temperature(scenario)

    u₁₀ = _windspeed(scenario, 10.0)
    ρₐ = _atmosphere_density(scenario)
    Tₐ = _atmosphere_temperature(scenario)

    # Setting up the Britter-McQuaid curves
    britter_interps = [ ]
    concs = sort(collect(keys(Britter_McQuaid_correlations)), rev=true)

    for conc in concs
        αs, βs = Britter_McQuaid_correlations[conc]
        f = extrapolate(interpolate((αs,), βs, Gridded(Linear())), Line())
        push!(britter_interps, (c=conc, it=f))
    end

    # relative density
    g = 9.80616  # gravity, m/s^2
    gₒ = g * ((ρⱼ - ρₐ)/ ρₐ)

    # critical distance
    Vr = Q/ρⱼ # volumetric release rate
    D = √(Vr/u₁₀)

    # temperature correction
    T′ = Tₐ/Tᵣ

    # correlation parameter
    α = 0.2*log10( gₒ^2 * Vr * u₁₀^-5 )

    # setting up correlation for constant α
    # calculates the points for the linear interpolation
    concs = [ elem.c for elem in britter_interps ]
    βs = [ elem.it(α) for elem in britter_interps ]

    # linear interpolation between the short distance correlation
    # and the main correlation
    βs = [ log10(30); βs]
    concs = [ 306*30^-2 / (1+ 306*30^-2); concs ]

    # linear interpolation, extrapolates past the end with a straight line
    itp = extrapolate(interpolate((βs,), concs, Gridded(Linear())), Line())

    # plume dimension parameters
    lb  = Vr*gₒ/(u₁₀^3)

    return BritterMcQuaidPlumeSolution(
        scenario, #scenario::Scenario
        :brittermcquaid, #model::Symbol
        ρⱼ,    # jet_density
        T′,    # temperature_correction
        D,     # critical length, m
        lb,    # length parameter, m
        itp    # interpolation::Extrapolation
    )
end

function (b::BritterMcQuaidPlumeSolution)(x, y, z)
    D = b.D
    lb = b.lb
    Lu  = 0.5*D + 2*lb
    Lh0 = D + 8*lb
    # domain check
    if x < -Lu
        return 0.0
    end

    # plume dimensions
    if x < 0
        # connecting the upwind extent to the crosswind extent at x=0
        # using a similar curve as for x>0
        x′ = x + Lu
        A = (Lh0^3)/(Lu^2)
        Lh = ∛(A*x′^2)
        Lv = (D^2)/Lh0
    else
        Lh = Lh0 + 2.5∛(lb*x^2)
        Lv = (D^2)/Lh0
    end

    # within the extent of the plume
    if (y < -Lh) || (y > Lh) || (z < 0) || (z > Lv)
        return 0.0
    end

    x′ = x/D
    if x′ < 30
        # use short distance correlation
        c′ = x′ > 0 ? 306*x′^-2 / (1+ 306*x′^-2) : 1.0
    else
        # use linear interpolation
        β = log10(x′)
        c′ = b.itp(β)
    end
    c = ( b.ρⱼ*c′*b.T′) / (1 - c′ + c′*b.T′)
    return c
end
