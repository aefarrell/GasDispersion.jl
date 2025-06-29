_mixing_height(s::Scenario) = _mixing_height(s.atmosphere)

_mixing_height(a::SimpleAtmosphere) = Inf # just for a test