# original SLAB wind profile function
function _slab_uafn(z,z0,ala0,zl,hmx,zt,cu1,cu2)
    #c  calculates the normalized ambient wind speed at height z

    if z < zt
        uatm = cu1*z + cu2*z*z
    else
        if ala0 < 0
            phmi = 1/sqrt(sqrt(1 - 16*zl*ala0))
            gus = -8*ala0/(1 - phmi)
            xu = sqrt(1 + gus*z)
            xu0 = sqrt(1 + gus*z0)
            uatm = log(z/z0) - phmi*(z-z0)/hmx - 2*(1 - phmi)*
                  (log((1 + xu)/(1 + xu0)) - (sqrt(xu) - sqrt(xu0))/
                  (gus*hmx))
        else
            uatm = log(z/z0) - (z-z0)/hmx + 5*ala0 * zl * ((1+zl/hmx)*
                  log((z+zl)/(z0+zl)) - (z-z0)/hmx)
        end
    end
    return uatm
end

_slab_uafn(z,wp::SLAB_Wind_Profile) = _slab_uafn(z,wp.z0,wp.ala0,wp.zl,
                                                 wp.hmx,wp.zt,wp.cu1,
                                                 wp.cu2)