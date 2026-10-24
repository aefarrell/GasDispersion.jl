# SLAB uses the Antoine equation for saturation pressures
# it also uses the derivative of the saturation pressure curve
# for various numerical routines

struct SLAB_Antoine{T <: Number}
    A::T
    B::T
    C::T
end

function _slab_psat(T, antoine::SLAB_Antoine)
    # calculates the saturation vapor pressure of water at temperature T
    # using the Antoine equation
    A = antoine.A
    B = antoine.B
    C = antoine.C
    return exp(A - B/(T + C))
end

function _slab_∂psat_∂T(T, antoine::SLAB_Antoine)
    # calculates the derivative of the saturation vapor pressure of water
    # with respect to temperature at temperature T using the Antoine equation
    A = antoine.A
    B = antoine.B
    C = antoine.C
    return _slab_psat(T, antoine) * B / (T + C)^2
end