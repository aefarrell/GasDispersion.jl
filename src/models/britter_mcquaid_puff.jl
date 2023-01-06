struct BritterMcQuaidPuff <: PuffModel end

"""
    puff(scenario::Scenario, BritterMcQuaidPuff)

Generates a Britter-McQuaid dispersion model on the given scenario and returns a
callable struct giving the concentration of the form
c(x, y, z, t)

"""
function puff(scenario::Scenario, ::Type{BritterMcQuaidPuff})
    # TODO
    error("Britter-McQuaid puff dispersion model is not currently implemented")
end
