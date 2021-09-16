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
    britter_plume_factory(scenario)

Generates a Britter-McQuaid dispersion model on the given scenario and returns a
function giving the centerline concentration of the form
c(x, [y, z, t])

Currently only implements the max concentration at a downwind distance x, the
other coordinates are ignored.

"""
function britter_plume_factory(scenario)
    Q = scenario.mass_emission_rate
    h = scenario.release_height
    ρⱼ = scenario.jet_density
    Tᵣ = scenario.release_temperature

    u = scenario.windspeed
    ρₐ = scenario.ambient_density
    Tₐ = scenario.ambient_temperature
    class = scenario.pasquill_gifford

    # Setting up the Britter-McQuaid curves
    britter_interps = [ ]
    concs = sort(collect(keys(Britter_McQuaid_correlations)), rev=true)

    for conc in concs
        αs, βs = Britter_McQuaid_correlations[conc]
        f = LinearInterpolation(αs, βs, extrapolation_bc=Line())
        push!(britter_interps, (c=conc, it=f))
    end


    # determining the windspeed at 10m
    wind = windspeed(u, h, class)
    u₁₀ = wind(10)

    # relative density
    g = 9.80616  # gravity, m/s^2
    gₒ = g * ((ρⱼ - ρₐ)/ ρₐ)

    # critical distance
    Vr = Q/ρⱼ # volumetric release rate
    D = √(Vr/u₁₀)

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
    it_new = LinearInterpolation(βs, concs, extrapolation_bc=Line())

    function britter_plume(x, y=0, z=0)
        T′ = Tₐ/Tᵣ
        x′ = x/D
        if x′ < 30
            # use short distance correlation
            c′ = x′ > 0 ? 306*x′^-2 / (1+ 306*x′^-2) : 1.0
        else
            # use linear interpolation
            β = log10(x′)
            c′ = it_new(β)
        end
        c = ( ρⱼ*c′*T′) / (1 - c′ + c′*T′)
        return c
    end

    return britter_plume
end
