
struct SimpleJet <: PlumeModel
    release_angle::Number
    k2::Number
    k3::Number
end
SimpleJet(; release_angle=0, k2=6, k3=5) = SimpleJet(release_angle,k2,k3)

struct SimpleJetSolution <: Plume
    scenario::Scenario
    model::Symbol
    release_angle::Number
    initial_concentration::Number
    rotation_matrix::Matrix
    k2::Number
    k3::Number
end

@doc doc"""
    plume(scenario::Scenario, SimpleJet(release_angle=0.0, k2=6, k3=5)

Generates a simple turbulent jet dispersion model on the given scenario and
returns a callable struct giving the concentration of the form
c(x, y, z[, t])

```math
c\left(x,y,z\right) = k_2 c_0 \left( d \over z \right) \sqrt{ \rho_j \over \rho_a }
\exp \left( - \left( k_2 { y \over x } \right)^2 \right)
\left[ \exp \left( - \left( k_2 { (z-h) \over x }\right)^2 \right)
+ \exp \left( - \left( k_3 { (z+h) \over x }\right)^2 \right) \right]
```

Assumes a circular jet with diameter equal to the jet diameter. Uses a gaussian
concentration profile.

`release_angle` controls the angle at which the jet is released, in radians
`k2` and `k3` are the parameters of the turbulent jet model.
"""
function plume(scenario::Scenario, model::SimpleJet)

    required_params = [:mass_emission_rate, :jet_diameter, :jet_velocity, :jet_density,
                       :ambient_density, :release_height]
    if all(key -> !(ismissing(getproperty(scenario,key))), required_params)
        # Rotation matrix
        θ = -1*model.release_angle
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

    return SimpleJetSolution(
    scenario,      #scenario::Scenario
    :simple_jet,   #model::Symbol
    model.release_angle, #release_angle::Number
    c₀,            #initial_concentration::Number
    R,             #rotation_matrix::Matrix
    model.k2,
    model.k3
    )

end

function (j::SimpleJetSolution)(x, y, z, t=0)
    d  = j.scenario.jet_diameter
    uⱼ = j.scenario.jet_velocity
    ρⱼ = j.scenario.jet_density
    ρₐ = j.scenario.ambient_density
    h  = j.scenario.release_height
    c₀ = j.initial_concentration
    Rθ = j.rotation_matrix
    k₂ = j.k2
    k₃ = j.k3

    # rotated jet
    x′, y′, z′ = Rθ*[x y z-h]'
    if x′ ≤ 0
        c = 0.0
    else
        ξ² = (y′^2 + z′^2)/(x′^2)
        c  = c₀*k₂*√(ρⱼ/ρₐ)*(d/x′)*exp(-k₃^2*ξ²)
    end

    # reflected jet (method of images)
    xᵣ, yᵣ, zᵣ = Rθ*[x y -z-h]'
    if xᵣ ≤ 0
        cᵣ = 0.0
    else
        ξᵣ² = (yᵣ^2 + zᵣ^2)/(xᵣ^2)
        cᵣ  = c₀*k₂*√(ρⱼ/ρₐ)*(d/xᵣ)*exp(-k₃^2*ξᵣ²)
    end

    return min(c+cᵣ, c₀)

end
