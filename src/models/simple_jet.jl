
struct SimpleJet <: PlumeModel
    release_angle::Number
    k2::Number
    k3::Number
end
SimpleJet(; release_angle=0, k2=6, k3=5) = SimpleJet(release_angle,k2,k3)

struct SimpleJetSolution <: Plume
    scenario::Scenario
    model::Symbol
    diameter::Number
    height::Number
    angle::Number
    initial_concentration::Number
    density_correction::Number
    k2::Number
    k3::Number
end

@doc doc"""
    plume(scenario::Scenario, SimpleJet(kwargs...))

Generates a simple turbulent jet dispersion model for the given scenario and
returns a callable giving the concentration with the form `c(x, y, z[, t])`

Turbulent jet model is per Long, V.D., "Estimation of the Extent of Hazard Areas
Around a Vent", *Chem. Process Hazard*, II, 6, 1963 ↩

```math
c\left(x,y,z\right) = k_2 c_0 \left( d \over z \right) \sqrt{ \rho_j \over \rho_a }
\exp \left( - \left( k_3 { r \over z } \right)^2 \right)
```

where *r* is the radial distance from the jet centerline. Assumes a circular jet
with diameter equal to the jet diameter.Ground-reflection is included by method
of images.

# Arguments
- `release_angle::Number=0`: the angle at which the jet is released, in radians
- `k2::Number=6` parameter of the model, default value is recommended by Long
- `k3::Number=5` parameter of the model, default value is recommended by Long
"""
function plume(scenario::Scenario, model::SimpleJet)
    # Density correction
    ρj = _release_density(scenario)
    ρa = _atmosphere_density(scenario)
    kd = √(ρj/ρa)

    # Initial concentration
    m = _mass_rate(scenario)
    d = _release_diameter(scenario)
    Q = _release_flowrate(scenario)
    c0 = m/Q

    return SimpleJetSolution(
    scenario,      #scenario::Scenario
    :simple_jet,   #model::Symbol
    d,  # diameter
    _release_height(scenario),  # height
    -1*model.release_angle, # release angle
    c0, # concentration
    kd, # density_correction
    model.k2,
    model.k3
    )

end

function (j::SimpleJetSolution)(x, y, z, t=0)
    d  = j.diameter
    h  = j.height
    θ  = j.angle
    c0 = j.initial_concentration
    kd = j.density_correction
    k₂ = j.k2
    k₃ = j.k3

    # rotated jet
    x′ = x*cos(θ)-(z-h)*sin(θ)
    y′ = y
    z′ = x*sin(θ)+(z-h)*cos(θ)
    if x′ ≤ 0
        c = 0.0
    else
        ξ² = (y′^2 + z′^2)/(x′^2)
        c  = c0*k₂*kd*(d/x′)*exp(-k₃^2*ξ²)
    end

    # reflected jet (method of images)
    xᵣ = x*cos(θ)-(-z-h)*sin(θ)
    yᵣ = y
    zᵣ = x*sin(θ)+(-z-h)*cos(θ)
    if xᵣ ≤ 0
        cᵣ = 0.0
    else
        ξᵣ² = (yᵣ^2 + zᵣ^2)/(xᵣ^2)
        cᵣ  = c0*k₂*kd*(d/xᵣ)*exp(-k₃^2*ξᵣ²)
    end

    return min(c+cᵣ, c0)

end
