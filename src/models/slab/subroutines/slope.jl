function _slab_sub_slope!(f,p::SLAB_Params,rho,x,h,v,w,b,bb,vg,u,wc,cm,ft,fu,fv,fw,bse)
    # unpacking parameters
    qs   = p.spl.qs
    rhoa = p.met.rhoa
  
    #subroutine slope
    #c  calculates the gradients used in the steady state plume conservation
    #c   equations

    #c  mass equation

    if x >= bse
        rqs = 0.0
    else
        rqs = 0.25*qs/bse
    end

    vhwb = √(3)*(v*h + w*bb)
    f[1] = rhoa*vhwb + rqs

    #c  half-width equation
    f[2] = (√(3)*(rhoa/rho)*v + vg)/u

    #c  b parameter equation
    f[3] = (vg*b)/(u*bb)

    #c  crosswind velocity equation
    f[4] = alfgv*grav*(rho-rhoa)*h*h
    f[10] = fv

    #c  temperature flux equation
    f[5] = ft

    #c  downwind velocity flux equation
    f[6] = fu

    #c  vertical velocity equation
    f[7] = -grav*(rho-rhoa)*bb*h
    f[11] = fw

    #c  center height equation
    f[8] = wc/u

    #c  released mass in cloud
    f[9] = 2*rho*bb*h*cm

    #end slope
end