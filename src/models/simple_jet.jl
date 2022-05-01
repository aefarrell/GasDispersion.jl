struct SimpleJet <: PlumeModel
    scenario::Scenario
    model::Symbol
    jetmodel::Symbol
    release_angle::Number
    initial_concentration::Number
    rotation_matrix::Matrix
    k₁::Number
    k₂::Number
end

"""
    simple_jet_factory(scenario; release_angle=0.0, jet_model=:gaussian, k₁=6, k₂=5)

Generates a simple turbulent jet dispersion model on the given scenario and
returns a callable struct giving the concentration of the form
c(x, y, z[, t])

Assumes a circular jet with diameter equal to the jet diameter.

`release_angle` controls the angle at which the jet is released, in radians
`jet_model` controls the turbulent jet model, by default a gaussian
`k₁` and `k₂` are the parameters of the turbulent jet model.
"""
function simple_jet_factory(scenario::Scenario; release_angle=0.0, jetmodel=:gaussian, k₁=6, k₂=5)

    required_params = [:mass_emission_rate, :jet_diameter, :jet_velocity, :jet_density,
                       :ambient_density, :release_height]
    if all(key -> !(ismissing(getproperty(scenario,key))), required_params)
        # Rotation matrix
        θ = -1*release_angle
        R = [cos(θ) 0 -sin(θ); 0 1 0; sin(θ) 0 cos(θ)]

        # Initial concentration
        m = scenario.mass_emission_rate
        d = scenario.jet_diameter
        u = scenario.jet_velocity
        Q = u*(π/4)*d^2 # volumetric flow rate
        c₀ = m/Q

    else
        missing_params = [ String(i) for i in filter(key -> ismissing(getproperty(scenario,key)), required_params)]
        error_string = "These parameters cannot be missing: " * join(missing_params, ", ")
        e = MissingException(error_string)
        throw(e)
    end

    return SimpleJet(
    scenario,      #scenario::Scenario
    :gaussian,     #model::Symbol
    jetmodel,      #jetmodel::Symbol
    release_angle, #release_angle::Number
    c₀,            #initial_concentration::Number
    R,             #rotation_matrix::Matrix
    k₁,
    k₂
    )

end

function (j::SimpleJet)(x, y, z, t=0)
    d  = j.scenario.jet_diameter
    uⱼ = j.scenario.jet_velocity
    ρⱼ = j.scenario.jet_density
    ρₐ = j.scenario.ambient_density
    h  = j.scenario.release_height
    c₀ = j.initial_concentration
    Rθ = j.rotation_matrix
    k₁ = j.k₁
    k₂ = j.k₂

    # rotated jet
    x′, y′, z′ = Rθ*[x y z-h]'
    if x′ ≤ 0
        c = 0.0
    else
        ξ² = (y′^2 + z′^2)/(x′^2)
        c  = c₀*k₁*√(ρⱼ/ρₐ)*(d/x′)*exp(-k₂^2*ξ²)
    end

    # reflected jet (method of images)
    xᵣ, yᵣ, zᵣ = Rθ*[x y -z-h]'
    if xᵣ ≤ 0
        cᵣ = 0.0
    else
        ξᵣ² = (yᵣ^2 + zᵣ^2)/(xᵣ^2)
        cᵣ  = c₀*k₁*√(ρⱼ/ρₐ)*(d/xᵣ)*exp(-k₂^2*ξᵣ²)
    end

    return min(c+cᵣ, c₀)

end
