function _slab_sub_entran(p::SLAB_Params,idpf,xn,timn,wss,xstr,ubs20,u,ug,vg,uab,rho,
                          zc,t,h,htp,bb,bbx,wc,_cp)
    # unpacking parameters
    idspl = p.spl.idspl
    rhoa = p.met.rhoa
    uastr = p.met.uastr
    ta = p.met.ta
    ala = p.met.ala
    z0 = p.met.z0
    ala0 = p.wps.ala0
    zl = p.wps.zl
    hmx = p.wps.hmx
    tgon = p.ocs.tgon
    bse = p.ocs.bse
    urf = p.ocs.urf
    rcf = p.ocs.rcf
    afa = p.ocs.afa
    stb = p.ocs.stb
    phimi = p.ocs.phimi
    phgam = p.ocs.phgam

    #subroutine entran

    #c  vertical entrainment
    
    ubar = .5*ug
    urho = ubar*rhoa/rho
    vbar = .5*vg
    vrho = vbar*rhoa/rho
    delu = (uab - u) * (rhoa / rho)
    cf = uastr/uab
    umgs = cf * u
    ugmgs = cf * ubar
    vmgs = cf * vbar

    uhs2 = cf1*(delu*delu + vrho*vrho)
    
    if htp > h
        ubsi2 = uastr^2 + uhs2
    else
        ubsi2 = umgs*umgs + vmgs*vmgs + cws*wss*uab
    end

    ubs2 = ubsi2
    if idspl == 1
        if xn > bse
            ubs2 = ubsi2 + (ubs20-ubsi2)*exp((bse-xn)/xstr)
        end
    end

    if idspl == 4
        qs = p.spl.qs
        qtcs = p.spl.qtcs
        tsd = p.spl.tsd
        if (qs*timn) > qtcs
            ubs2 = ubsi2 + (ubs20-ubsi2)*exp((tsd-timn)*u/xstr)
        end
    end

    ubs = sqrt(ubs2)
    
    vh = cf*ubs
    tg = tgon*ta + (1 - tgon)*t

    if htp > h
        uts2 = 0.0
    else
        if tg > t
            tgt = 0.5*(tg + t)
            uts2 = ((cth*grav*(tg-t)*vh*h/tgt)^tutrd)
        else
            uts2 = 0.0
        end
    end
    
    
    #2045 
    ustr2 = ubs2 + uhs2 + uts2
    ustr = sqrt(ustr2)
    
    galc = cri*grav*(rho - rhoa)/rho
    if galc < 0
        galc = 0.0
    end

    alah = ala0/(1 + htp/zl)
    al = (alah*uastr^2 + galc)/ustr2

    if al < 0
        phiit = 1/sqrt(1 - 16*htp*al)
    else
        phiit = 1 + 5*htp*al
    end
    
    fpr = 1 - htp/hmx
    fprb = 0.0
    phiitb = 1.0
    htst = 1.01*h

    if htp > htst
        hbt = htp - h
        fprb = 1 - hbt/hmx

        if al < 0
            phiitb = 1/sqrt(1 - 16*hbt*al)
        else
            phiitb = 1 + 5*hbt*al
        end
    end

    uah = _slab_uafn(htp,p.wps)
    w = (urf/uah)*aa*xk*ustr*(fpr/phiit + fprb/phiitb)
    
    #c  horizontal entrainment
    
    #al = ala this seems pointless?
    if ala < 0
        stby = 1 - rcf*10*ala
    else
        stby = 1/(1 + rcf*10*ala)
    end

    atot = afa*stby
    rab = 0.5*sigb/atot
    va = atot*uab/(1 + rab*bb/sr3)
    vjp = aa*xk*delu
    v = sqrt(va*va + cf1*vjp*vjp)
    
    if zc > (0.5*h)
        sig = 0.5*h/sr3
    else
        sig = (h-zc)/sr3
    end

    zrf = zc + 0.5*sig
    frzr = zrf/(zrf+z0)
    fgr = 1 - zrf/hmx

    if stb < 0
        phmr = phimi + (1 - phimi)/sqrt(1 + phgam*zrf)
    else
        phmr = 1 + (5*ala0*zrf)/(1 + zrf/zl)
    end
    
    vxs = 0.6*uastr*frzr*phmr*fgr/xk
    vax = atot*uab/(1 + rab*bbx/sr3)
    vx = sqrt(vax*vax+vxs*vxs)
    
    #c  velocity and temperature fluxes
    
    # I don't think this does anything?
    # if ala < 0
    #     phiit = 1/sqrt(1 - 16*htp*ala)
    # else
    #     phiit = 1 + 5*htp*ala
    # end

    umhs2 = cf1*abs(delu)*delu
    ugmhs2 = cf1*abs(urho)*urho
    vmhs2 = cf1*abs(vrho)*vrho
    wmhs2 = cf1*abs(wc)*wc
    rhoef = rho

    if idpf == 2
        rhoef = rho*bbx
    end

    if htp > h
        fu = rhoef*(bb+bb+h)*umhs2
        fug = 0.0
        fv = 0.0
        fw = -rhoef*h*wmhs2
        ft = 0.0
    else
        fu = -rhoef*(bb*(umgs*umgs - uastr^2)-(bb + h)*umhs2)
        fug = -rhoef*bb*(ugmgs*ugmgs + ugmhs2)
        fv = -rhoef*bb*(vmgs*vmgs + vmhs2)
        fw = 0.0
        ft = rhoef*bb*vh*_cp*(tg - t)
    end

    #end entran
    
    return w,v,vx,ubs2,fug,ft,fu,fv,fw
end