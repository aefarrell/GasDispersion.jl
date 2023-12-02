function _slab_sub_evalpf(p::SLAB_Params,sru0,zc,h0,rho,rho0,bb0,b0,bbx0,bx0,r,r0,
                          bbv,bbv0,bv,bv0,bbvx,bbvx0,bvx,bvx0,u0,uab0,cm,vg0,ug0,
                          wc0,htp,htp0,sfu,sfx,sfy,sfz,g,gw,gx)
    # unpacking parameters
    uastr = p.met.uastr
    rhoa = p.met.rhoa
    hmx = p.wps.hmx

    #subroutine evalpf

    #c  calculates the transient puff velocities and height
    
    #c  coupled velocity, height and average ambient velocity calculation

    u1 = (uastr/xk)*_slab_uafn(htp0,p.wps)
    gey = bbx0*h0/(h0*(bbx0+bb0) + 2*bb0*bbx0)
    gex = (bb0/bbx0)*gey
    rey = (rho0/rho)^gey
    rex = (rho0/rho)^gex
    bb = (bb0+bbv-bbv0)*rey
    b = (b0+bv-bv0)*rey
    bbx = (bbx0+bbvx-bbvx0)*rex
    bx = (bx0+bvx-bvx0)*rex

    beta = sqrt(bb*bb-b*b)/√(3)
    betax = sqrt(bbx*bbx-bx*bx)/√(3)

    h = r/(rho*bbx*bb)
    if h > hmx
        h = hmx
    end

    hhf = 0.5*h
    zcmx = hmx-hhf

    if zc > zcmx
        zc = zcmx
    end

    if zc > hhf
        htp = zc+hhf
    else
        htp = h
    end

    htph = .5*(htp0+htp)
    u12 = (uastr/xk)*_slab_uafn(htph,p.wps)
    u2 = (uastr/xk)*_slab_uafn(htp,p.wps)
    uab = (htp0*uab0+(u1 + 4*u12+u2)*(htp-htp0)/6.)/htp

    if sfu == 0.0
        sfue = 0.0
    else
        sfue = -r0*(u0-uab0)*abs(sfu)/(abs(sfu) + r0*abs(u0-uab0))
    end

    u = (1 - cm)*uab + (sru0+sfue)/r

    #c  calculate crosswind velocities

    if htp > h
        vg = 0.0
        ug = 0.0
        if wc0 == 0
            ew = 1.0
        else
            ew = (r0*abs(wc0))/(r0*abs(wc0) + abs(sfz))
        end
        wc = ew*(gw + r0*wc0)/r
    else
        if rho > rhoa
            if vg0 == 0.
                ev = 1.0
                rbxy = bb0*bbx0/(bb0*bb0+bbx0*bbx0)
                if h0 > 0.
                    vg0 = -2*rbxy*bbx0*wc0/h0
                end
            else
                ev = (r0*vg0)/(r0*vg0 + abs(sfy))
            end
            vg = ev*(g + r0*vg0)/r
            if ug0 == 0
                eu = 1.0
                rbxy = bb0*bbx0/(bb0*bb0+bbx0*bbx0)
                if h0 > 0
                    ug0 = -2*rbxy*bb0*wc0/h0
                end
            else
                eu = (r0*ug0)/(r0*ug0 + abs(sfx))
            end
            ug = eu*(gx + r0*ug0)/r
            wc = -(vg/bb+ug/bbx)*zc
        else
            if wc0 == 0
                ew = 1.0
            else
                ew = (r0*abs(wc0))/(r0*abs(wc0) + abs(sfz))
            end
            wc = ew*(gw + r0*wc0)/r
            vg = 0.0
            ug = 0.0
        end
    end

    return bb,b,bbx,bx,beta,betax,h,zc,htp,u,uab,vg,vg0,ug,ug0,wc,wc0
end