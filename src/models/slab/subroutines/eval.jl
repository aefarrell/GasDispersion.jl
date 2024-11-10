function _slab_sub_eval(p::SLAB_Params,xn,alfg,sru0,zc,h0,u0,uab0,b0,bb0,r,r0,bv,
                         bv0,bbv,bbv0,rho,rho0,vg0,wc0,cm,htp0,b,htp,uab,beta,vg,
                         wc,h,u,bb,sfu,sfz,sfy,g,gw,bse)
    # unpacking parameters
    uastr = p.met.uastr
    rhoa = p.met.rhoa
    hmx = p.wps.hmx

    #subroutine eval
    #c  calculates the steady state plume velocities and height
    #c  coupled velocity, height and average ambient velocity calculation

    # these variables are used before being assigned
    utc3 = 0.0
    utm3 = 0.0

    #c   calculate downwind velocity
    u1 = (uastr/xk)*_slab_uafn(htp0,p.wps)
    u = u0
    bb = bb0 + bbv - bbv0
    h = r/(rho*u*bb)
    rex = (rho0/rho)^(h/(h+bb+bb))

    #do 1050 i=1,11
    for i = 1:11
        gex = h/(h+bb+bb)
        if u0 > uab0
            rex0 = rho0*u0/(rho*u)
        else
            rex0 = rho0/rho
        end

        rex = sqrt(rex*(rex0^gex))
        bb = (bb0 + bbv - bbv0)*rex
        b = (b0 + bv - bv0)*rex
        h = r/(rho*u*bb)

        if h > (bb+bb)
            bb = sqrt(r/(2 *rho*u))
            h = bb+bb
            rex = sqrt(rex0)
        end

        hhf = 0.5*h
        if zc > hhf
            htp = zc + hhf
        else
            htp = h
        end

        htph = 0.5*(htp0 + htp)
        u12 = (uastr/xk)*_slab_uafn(htph,p.wps)
        u2 = (uastr/xk)*_slab_uafn(htp,p.wps)
        uab = (htp0*uab0+(u1+4*u12+u2)*(htp-htp0)/6)/htp

        if uab > u2
            uab = u2
        end

        if sfu == 0
            sfue = 0.0
        else
            sfue = -r0*(u0-uab0)*abs(sfu)/(abs(sfu) + r0*abs(u0-uab0))
        end

        ut0 = uab*(1-cm) + (sru0 + sfue)/r
        utg3 = alfg * grav * (rho-rhoa)*r/(2*bb*rho*rho)

        if utg3 <= 0
            ut2 = ut0
            #go to 1040
        else
            #1015 
            utm3 = (2/27)*ut0^3
            utc3 = utg3 - utm3

            if utc3 > utm3
                ut2 = (2/3)*ut0
                ut2 = min(ut2,u0)
                #go to 1040
            else
                #1020 
                angl = acos(-utc3/utm3)
                ut2 = (ut0 + 2*ut0*cos(angl/3))/3
            end
        end

        #1040
        ut2 = (u*ut2^2)^.333333
        dluk = abs(u-ut2)/sqrt(u*ut2)
        u = ut2
        h = r/(rho*u*bb)

        if (dluk < 0.001) 
            #go to 1052
            break
        end
    #1050 continue
    end

    #1052 continue

    if (h < (zc+zc)) && ((alfg*rho) > (alfg*rhoa))
        h = zc + zc
        u = r/(rho*bb*h)
    end

    if (h > hmx) 
        h = hmx
    end

    hhf = .5*h
    zcmx = hmx - hhf

    if zc > zcmx
        zc = zcmx
    end

    if zc > hhf
        htp = zc + hhf
    else
        htp = h
    end

    htph = .5*(htp0+htp)
    u12 = (uastr/xk)*_slab_uafn(htph,p.wps)
    u2 = (uastr/xk)*_slab_uafn(htp,p.wps)
    uab = (htp0*uab0+(u1+4*u12+u2)*(htp-htp0)/6)/htp

    if uab > u2
        uab = u2
    end

    # dbbx = 0.0 assigned but never used
    
    # only used for evaporating pool releases
    # if utc3 > utm3
    #     if keval == 4
    #         # utst = 3.0 assigned but never used
    #         ru3 = (27.0/4.0)*utg3/(ut0*ut0*ut0)
    #         bb = ru3*bb
    #         b = ru3*b
    #         # dbbx = bb*(ru3 - 1)/ru3 assigned but never used
    #     end
    # end

    beta = sqrt(bb*bb-b*b)/√(3)

    #c  calculate crosswind velocities

    if xn > bse
        if htp > h
            vg = 0.0

            if wc0 == 0
                ew = 1.0
            else
                ew = r0*abs(wc0)/(r0*abs(wc0) + abs(sfz))
            end

            wc = ew*(gw + r0*wc0)/r
        else
            if rho > rhoa

                if vg0 == 0
                    ev = 1.0
                    vg0 = -2*wc0*bb0/h0
                else
                    ev = r0*vg0/(r0*vg0 + abs(sfy))
                end

                vg = ev*(g + r0*vg0)/r
                wc = -vg*zc/bb
            else

                if wc0 == 0
                    ew = 1.0
                else
                    ew = r0*abs(wc0)/(r0*abs(wc0) + abs(sfz))
                end

                wc = ew*(gw + r0*wc0)/r
                vg = 0.0
            end
        end

    else
        if rho > rhoa

            if vg0 == 0
                ev = 1.0
            else
                ev = r0*vg0/(r0*vg0 + abs(sfy))
            end
            
            vg = ev*(g + r0*vg0)/r
        else
            vg = 0.0
        end
        wc = 0.0
    end
    #end eval
    return u,uab,b,bb,beta,h,zc,vg,vg0,wc,htp
end