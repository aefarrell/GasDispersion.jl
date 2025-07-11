# defining type for dispatch
struct SimpleJet <: PlumeModel end

struct SimpleJetSolution{F<:Number} <: Plume
    scenario::Scenario
    model::Symbol
    diameter::F
    height::F
    angle::F
    initial_concentration::F
    density_correction::F
    k2::F
    k3::F
end
SimpleJetSolution(s,m,d,h,θ,c0,kd,k2,k3) = SimpleJetSolution(s,m,promote(d,h,θ,c0,kd,k2,k3)...)

_release_angle(::HorizontalJet) = 0.0
_release_angle(::VerticalJet) = π/2

# for reverse compatibility
function plume(scenario::Scenario, ::Type{SimpleJet}, eqs=DefaultSet; kwargs...) 
    @warn "plume(scenario, SimpleJet, eqs) is deprecated, use plume(scenario, SimpleJet(), eqs) instead."    
    return plume(scenario, SimpleJet(), eqs; kwargs...)
end

@doc doc"""
    plume(::Scenario, ::SimpleJet; kwargs...)

Returns the solution to a simple turbulent jet dispersion model for the given
scenario.

```math
c\left(x,y,z\right) = k_2 c_0 \left( d \over z \right) \sqrt{ \rho_j \over \rho_a }
\exp \left( - \left( k_3 { r \over z } \right)^2 \right)
```

where *r* is the radial distance from the jet centerline. Assumes a circular jet
with diameter equal to the jet diameter. Ground-reflection is included by method
of images.

# References
+ Long, V.D. 1963. "Estimation of the Extent of Hazard Areas Round a Vent." *Chem. Process Hazard*. II:6

# Arguments
- `release_angle::Number=nothing`: the angle at which the jet is released, in radians, if nothing defaults to 0 for horizontal jets and π/2 for vertical jets
- `k2::Number=6` parameter of the model, default value is recommended by Long
- `k3::Number=5` parameter of the model, default value is recommended by Long

"""
function plume(scenario::Scenario, ::SimpleJet, eqs::EquationSet=DefaultSet; release_angle=nothing, k2=6.0, k3=5.0)
    # release angle
    if isnothing(release_angle)
        release_angle = _release_angle(scenario.release)
    end
    # Density correction
    ρj = _release_density(scenario)
    ρa = _atmosphere_density(scenario)
    kd = √(ρj/ρa)

    # Initial concentration
    m = _mass_rate(scenario)
    d = _release_diameter(scenario)
    Qt = _release_flowrate(scenario)
    Qi = m/ρj
    c0 = min(Qi/Qt,1.0)

    return SimpleJetSolution(
    scenario,      #scenario::Scenario
    :simple_jet,   #model::Symbol
    d,  # diameter
    _release_height(scenario),  # height
    -1*release_angle, # release angle
    c0, # concentration
    kd, # density_correction
    k2,
    k3
    )

end

function (j::SimpleJetSolution{F})(x, y, z, t=0) where {F}
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
