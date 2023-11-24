function _slab_sub_slopepf!(f,p::SLAB_Params,rqs,atau,x,u,ug,v,vg,vx,w,wc,b,bx,
                            bb,bbx,h,rho,fug,ft,fu,fv,fw)
    # unpacking parameters
    rhoa = p.met.rhoa

    #subroutine slopepf

    #c  calculates the gradients used in the transient puff convervation
    #c    equations

    #c  mass and downwind cloud center equations

    vhwb = sr3*((vx*bb+v*bbx)*h + w*bb*bbx)
    f[1] = rhoa*vhwb + rqs
    f[2] = u - atau*x

    #c  half-width equation

    f[3] = sr3*(rhoa/rho)*v + vg

    #c  b parameter equation

    f[4] = (vg*b)/bb

    #c  half-length equation
    f[5] = sr3*(rhoa/rho)*vx + ug

    #c  bx parameter equation

    f[6] = (ug*bx)/bbx

    #c  crosswind velocity equations

    fgrv = alfgv*grav*(rho-rhoa)*h*h

    f[7] = fgrv*bbx
    f[13] = fv
    f[8] = fgrv*bb
    f[14] = fug

    #c  temperature flux equation

    f[9] = ft

    #c  downwind velocity flux equation

    f[10] = fu

    #c  vertical velocity equation

    f[11] = -grav*(rho-rhoa)*bb*bbx*h
    f[15] = fw

    #c  center height equation

    f[12] = wc

end